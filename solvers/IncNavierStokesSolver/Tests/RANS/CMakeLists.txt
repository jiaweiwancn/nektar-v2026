# Small external build against an existing installation. It does not modify
# the shared Nektar++ build or installation while another experiment uses them.
cmake_minimum_required(VERSION 3.12)
project(RANSWallLayerPrototype CXX)
find_package(Nektar++ REQUIRED)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
get_filename_component(RANS_SOURCE_ROOT "${CMAKE_CURRENT_LIST_DIR}/../../../.." ABSOLUTE)
add_executable(RANSReducedModelChecks test_cls_epawf_reduced_model.cpp)
target_include_directories(RANSReducedModelChecks PRIVATE ${RANS_SOURCE_ROOT}/solvers)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(RANSReducedModelChecks PRIVATE -Wall -Wextra -Wpedantic -Wno-sign-compare -Werror)
endif()
add_executable(RANSKBoundQPChecks test_cls_epawf_k_bound_qp.cpp)
target_include_directories(RANSKBoundQPChecks PRIVATE ${RANS_SOURCE_ROOT}/solvers)
add_executable(RANSWallTimingChecks test_rans_wall_timing.cpp)
add_executable(RANSWallMonitorChecks test_rans_wall_monitor.cpp)
target_include_directories(RANSWallMonitorChecks PRIVATE ${RANS_SOURCE_ROOT}/solvers)
target_include_directories(RANSWallTimingChecks PRIVATE ${RANS_SOURCE_ROOT}/solvers)
enable_testing()
add_test(NAME RANSKBoundQPChecks COMMAND RANSKBoundQPChecks)
add_test(NAME RANSReducedModelChecks COMMAND RANSReducedModelChecks)
add_executable(WJMappedAlgebraChecks test_wj_mapped_algebra.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallFieldLayout.cpp)
target_include_directories(WJMappedAlgebraChecks PRIVATE ${RANS_SOURCE_ROOT}/solvers)
target_link_libraries(WJMappedAlgebraChecks ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
add_test(NAME WJMappedAlgebraChecks COMMAND WJMappedAlgebraChecks)
set_tests_properties(WJMappedAlgebraChecks PROPERTIES PROCESSORS 1 TIMEOUT 120
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
add_definitions(${NEKTAR++_DEFINITIONS} -DENABLE_NEKTAR_EXCEPTIONS)
include_directories(${RANS_SOURCE_ROOT}/library)
include_directories(SYSTEM ${NEKTAR++_INCLUDE_DIRS} ${NEKTAR++_TP_INCLUDE_DIRS})
add_library(RANSKernels STATIC
    ${RANS_SOURCE_ROOT}/library/SolverUtils/Turbulence/RANSModels.cpp
    ${RANS_SOURCE_ROOT}/library/SolverUtils/Turbulence/WJEARSM.cpp
    ${RANS_SOURCE_ROOT}/library/SolverUtils/Turbulence/RANSRobustVariant.cpp
    ${RANS_SOURCE_ROOT}/library/SolverUtils/Turbulence/RANSWallInnerStress.cpp
    ${RANS_SOURCE_ROOT}/library/SolverUtils/Turbulence/AnalyticalWallFunction.cpp
    ${RANS_SOURCE_ROOT}/library/SolverUtils/Turbulence/EquilibriumAWFProfiles.cpp
    ${RANS_SOURCE_ROOT}/library/SolverUtils/Turbulence/EquilibriumAWFTransition.cpp
    ${RANS_SOURCE_ROOT}/library/SolverUtils/Turbulence/RANSWallGeometry.cpp
    ${RANS_SOURCE_ROOT}/library/SolverUtils/Turbulence/RANSWallNormalMap.cpp)
target_compile_definitions(RANSKernels PRIVATE SOLVER_UTILS_EXPORTS)
add_executable(RANSAcceptedAWFProbe
    ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSAcceptedAWFProbe.cpp)
target_link_libraries(RANSAcceptedAWFProbe RANSKernels)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(RANSAcceptedAWFProbe PRIVATE -Wall -Wextra -Wpedantic -Wno-sign-compare -Werror)
endif()
add_test(NAME RANSAcceptedAWFKernel COMMAND RANSAcceptedAWFProbe
    --output ${CMAKE_CURRENT_BINARY_DIR}/rans-accepted-awf-kernel.json)
set_tests_properties(RANSAcceptedAWFKernel PROPERTIES PROCESSORS 1 TIMEOUT 180
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
add_executable(RANSWallLayerProbe
    ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSWallLayerProbe.cpp)
add_executable(RANSRobustVariantProbe
    ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSRobustVariantProbe.cpp)
target_link_libraries(RANSRobustVariantProbe RANSKernels)
add_executable(RANSWallInnerStressProbe
    ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSWallInnerStressProbe.cpp)
target_link_libraries(RANSWallInnerStressProbe RANSKernels)
add_executable(RANSKernelProbe
    ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSKernelProbe.cpp)
add_executable(RANSWallGeometryProbe
    ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSWallGeometryProbe.cpp)
add_executable(RANSCurvedNoteProbe
    ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSCurvedNoteProbe.cpp)
add_executable(RANSCurvedGeometryProbe
    ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSCurvedGeometryProbe.cpp)
add_executable(RANSObservationKernelProbe
    ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSObservationKernelProbe.cpp)
add_executable(RANSStressTangentAudit
    ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSStressTangentAudit.cpp)
add_executable(RANSConstitutiveProbe
    ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSConstitutiveProbe.cpp)
add_executable(RANSHeightReferenceChecks test_cls_epawf_height_reference.cpp)
target_link_libraries(RANSHeightReferenceChecks RANSKernels ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
add_executable(RANSProfileBranchChecks test_cls_epawf_profile_branches.cpp)
target_link_libraries(RANSProfileBranchChecks RANSKernels ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
add_test(NAME RANSProfileBranchChecks COMMAND RANSProfileBranchChecks)
add_executable(RANSEquilibriumProfilesProbe
    ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSEquilibriumProfilesProbe.cpp)
add_executable(RANSEquilibriumTransitionProbe
    ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSEquilibriumTransitionProbe.cpp)
target_link_libraries(RANSEquilibriumTransitionProbe RANSKernels)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(RANSEquilibriumTransitionProbe PRIVATE -Wall -Wextra -Wpedantic -Wno-sign-compare -Werror)
endif()
add_executable(RANSWallExpressionProbe ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSWallExpressionProbe.cpp)
target_link_libraries(RANSWallExpressionProbe ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(RANSWallExpressionProbe PRIVATE -Wall -Wextra -Wpedantic -Wno-sign-compare -Werror)
endif()
add_executable(RANSEquilibriumPatchProbe
    ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSEquilibriumPatchProbe.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallFieldLayout.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallPatch.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallInterface.cpp)
target_include_directories(RANSEquilibriumPatchProbe PRIVATE ${RANS_SOURCE_ROOT}/solvers)
target_link_libraries(RANSEquilibriumPatchProbe RANSKernels ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
add_executable(RANSEquilibriumTimeProbe
    ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSEquilibriumTimeProbe.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallFieldLayout.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallTimeIntegration.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSHistory.cpp)
target_include_directories(RANSEquilibriumTimeProbe PRIVATE ${RANS_SOURCE_ROOT}/solvers)
target_link_libraries(RANSEquilibriumTimeProbe RANSKernels ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
add_executable(RANSHistoryRejectionChecks
    test_cls_epawf_history_rejection.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallTimeIntegration.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSHistory.cpp)
target_include_directories(RANSHistoryRejectionChecks PRIVATE ${RANS_SOURCE_ROOT}/solvers)
target_link_libraries(RANSHistoryRejectionChecks ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
add_test(NAME RANSHistoryRejectionChecks COMMAND RANSHistoryRejectionChecks)
add_executable(RANSTransverseEvolutionProbe
    ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSTransverseEvolutionProbe.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallFieldLayout.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallTimeIntegration.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSHistory.cpp)
target_include_directories(RANSTransverseEvolutionProbe PRIVATE ${RANS_SOURCE_ROOT}/solvers)
target_link_libraries(RANSTransverseEvolutionProbe RANSKernels ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(RANSTransverseEvolutionProbe PRIVATE -Wall -Wextra -Wpedantic -Wno-sign-compare -Werror)
endif()
add_executable(RANSEquilibriumEvolutionProbe
    ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSEquilibriumEvolutionProbe.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallFieldLayout.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallPatch.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallGlobalField.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionBoundary.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionEvaluation.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionTransition.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionRobust.cpp)
target_include_directories(RANSEquilibriumEvolutionProbe PRIVATE ${RANS_SOURCE_ROOT}/solvers)
target_link_libraries(RANSEquilibriumEvolutionProbe RANSKernels ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
add_executable(RANSCurvedWallProbe
    ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSCurvedWallProbe.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallFieldLayout.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallPatch.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallGlobalField.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionBoundary.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionEvaluation.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionTransition.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionRobust.cpp)
target_include_directories(RANSCurvedWallProbe PRIVATE ${RANS_SOURCE_ROOT}/solvers)
target_link_libraries(RANSCurvedWallProbe RANSKernels ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
foreach(rans_native_target RANSEquilibriumEvolutionProbe RANSCurvedWallProbe)
    target_sources(${rans_native_target} PRIVATE
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionActions.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionNonpressure.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionResolution.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionFeasible.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPicard.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPimple.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionKBound.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionContact.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionSource.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionReduced.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionDiagnostics.cpp)
endforeach()
add_executable(RANSPressureActionProbe
    ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSPressureActionProbe.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallFieldLayout.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallPatch.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallGlobalField.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionBoundary.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionEvaluation.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionTransition.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionRobust.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionActions.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionNonpressure.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionResolution.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionFeasible.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPicard.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPimple.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionKBound.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionContact.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionSource.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionReduced.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionDiagnostics.cpp)
target_include_directories(RANSPressureActionProbe PRIVATE ${RANS_SOURCE_ROOT}/solvers)
target_link_libraries(RANSPressureActionProbe RANSKernels ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
foreach(rans_r3_probe RANSNonpressureActionProbe RANSWallResolutionProbe)
    add_executable(${rans_r3_probe}
        ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/${rans_r3_probe}.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallFieldLayout.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallPatch.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallGlobalField.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionBoundary.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionEvaluation.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionTransition.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionRobust.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionActions.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionNonpressure.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionResolution.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionFeasible.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPicard.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPimple.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionKBound.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionContact.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionSource.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionReduced.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionDiagnostics.cpp)
    target_include_directories(${rans_r3_probe} PRIVATE ${RANS_SOURCE_ROOT}/solvers)
    target_link_libraries(${rans_r3_probe} RANSKernels ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(${rans_r3_probe} PRIVATE -Wall -Wextra -Wpedantic -Wno-sign-compare -Werror)
    endif()
endforeach()
target_sources(RANSWallResolutionProbe PRIVATE
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallTimeIntegration.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSHistory.cpp)
add_executable(RANSDurableProgressProbe ${RANS_SOURCE_ROOT}/docs/plans/cls-epawf-revised-evidence/m8-r3/durable-probe.cpp)
target_include_directories(RANSDurableProgressProbe PRIVATE ${RANS_SOURCE_ROOT}/solvers ${RANS_SOURCE_ROOT}/library)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    foreach(rans_r2_target RANSWallInnerStressProbe RANSPressureActionProbe)
        target_compile_options(${rans_r2_target} PRIVATE -Wall -Wextra -Wpedantic -Wno-sign-compare -Werror)
    endforeach()
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    foreach(rans_target RANSRobustVariantProbe RANSKernels RANSWallLayerProbe RANSKernelProbe RANSWallGeometryProbe RANSCurvedGeometryProbe RANSCurvedNoteProbe RANSObservationKernelProbe RANSStressTangentAudit RANSConstitutiveProbe RANSEquilibriumProfilesProbe RANSEquilibriumPatchProbe RANSEquilibriumTimeProbe RANSEquilibriumEvolutionProbe RANSCurvedWallProbe)
        target_compile_options(${rans_target} PRIVATE
            -Wall -Wextra -Wpedantic -Wno-sign-compare -Werror)
    endforeach()
endif()
target_link_libraries(RANSWallLayerProbe
    RANSKernels ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
target_link_libraries(RANSKernelProbe RANSKernels)
target_link_libraries(RANSWallGeometryProbe RANSKernels)
target_link_libraries(RANSCurvedGeometryProbe RANSKernels)
target_link_libraries(RANSCurvedNoteProbe RANSKernels)
target_link_libraries(RANSObservationKernelProbe RANSKernels)
target_link_libraries(RANSStressTangentAudit RANSKernels)
target_link_libraries(RANSConstitutiveProbe RANSKernels)
target_link_libraries(RANSEquilibriumProfilesProbe RANSKernels ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
option(RANS_BUILD_SOLVER "Build the complete incompressible solver in isolation" OFF)
if(RANS_BUILD_SOLVER)
    set(RANS_PREBUILT_SPATIAL_LIBRARY "" CACHE FILEPATH
        "Optional isolated, source-compatible SpatialDomains library")
    if(RANS_PREBUILT_SPATIAL_LIBRARY)
        if(NOT EXISTS "${RANS_PREBUILT_SPATIAL_LIBRARY}")
            message(FATAL_ERROR "RANS_PREBUILT_SPATIAL_LIBRARY does not exist")
        endif()
        get_filename_component(rans_spatial_real "${RANS_PREBUILT_SPATIAL_LIBRARY}" REALPATH)
        get_filename_component(RANS_SPATIAL_RUNTIME_DIR "${rans_spatial_real}" DIRECTORY)
        add_library(RANSSpatialDomains SHARED IMPORTED)
        set_target_properties(RANSSpatialDomains PROPERTIES
            IMPORTED_LOCATION "${RANS_PREBUILT_SPATIAL_LIBRARY}")
    else()
    set(RANS_SPATIAL_RUNTIME_DIR "${CMAKE_CURRENT_BINARY_DIR}/geometry-lib")
    # Build the source SpatialDomains library in isolation too: terrain normals
    # depend on full-precision XML curve partitioning, which older installations
    # round. Use the production base list and installed dependency configuration.
    file(READ ${RANS_SOURCE_ROOT}/library/SpatialDomains/CMakeLists.txt rans_spatial_cmake)
    string(REGEX MATCH "SET\\(SPATIAL_DOMAINS_SOURCES[^)]*" rans_spatial_block "${rans_spatial_cmake}")
    string(REGEX MATCHALL "[A-Za-z0-9_/]+\\.cpp" rans_spatial_sources "${rans_spatial_block}")
    if(NEKTAR_USE_HDF5)
        list(APPEND rans_spatial_sources MeshGraphIOHDF5.cpp)
    endif()
    if(NEKTAR_USE_SCOTCH)
        list(APPEND rans_spatial_sources MeshPartitionScotch.cpp)
        if(NEKTAR_USE_MPI)
            list(APPEND rans_spatial_sources MeshPartitionPtScotch.cpp)
        endif()
    endif()
    if(NEKTAR_USE_METIS)
        list(APPEND rans_spatial_sources MeshPartitionMetis.cpp)
    endif()
    list(TRANSFORM rans_spatial_sources PREPEND "${RANS_SOURCE_ROOT}/library/SpatialDomains/")
    add_library(RANSSpatialDomains SHARED ${rans_spatial_sources})
    target_compile_definitions(RANSSpatialDomains PRIVATE SPATIAL_DOMAINS_EXPORTS)
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        # Match the production flags for existing Geometry sources.
        target_compile_options(RANSSpatialDomains PRIVATE -Wno-narrowing)
    endif()
    set_target_properties(RANSSpatialDomains PROPERTIES OUTPUT_NAME SpatialDomains
        VERSION ${NEKTAR++_VERSION} SOVERSION ${NEKTAR++_VERSION}
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/geometry-lib)
    endif()
    # Imported SpatialDomains exposes its partitioner dependencies. Retain those
    # while replacing only that library, rather than loading two copies.
    get_target_property(rans_spatial_dependencies SpatialDomains INTERFACE_LINK_LIBRARIES)
    if(RANS_PREBUILT_SPATIAL_LIBRARY)
        target_link_libraries(RANSSpatialDomains INTERFACE StdRegions ${rans_spatial_dependencies}
            ${NEKTAR++_TP_LIBRARIES})
    else()
        target_link_libraries(RANSSpatialDomains StdRegions ${rans_spatial_dependencies}
            ${NEKTAR++_TP_LIBRARIES})
    endif()
    set(RANS_SOLVER_ROOT ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver)
    # Use the production source list so registrations and legacy paths match
    # the normal executable. No shared installation is overwritten.
    file(READ ${RANS_SOLVER_ROOT}/CMakeLists.txt rans_solver_cmake)
    string(REGEX MATCH "LIBRARY_SOURCES[^)]*" rans_library_block "${rans_solver_cmake}")
    string(REGEX MATCHALL "\\./[A-Za-z0-9_/]+\\.cpp" rans_solver_sources "${rans_library_block}")
    set(rans_absolute_sources ${RANS_SOLVER_ROOT}/IncNavierStokesSolver.cpp)
    foreach(rans_source ${rans_solver_sources})
        list(APPEND rans_absolute_sources ${RANS_SOLVER_ROOT}/${rans_source})
    endforeach()
    add_executable(IncNavierStokesSolver ${rans_absolute_sources})
    target_include_directories(IncNavierStokesSolver PRIVATE ${RANS_SOURCE_ROOT}/solvers)
    target_link_libraries(IncNavierStokesSolver RANSKernels RANSSpatialDomains
        ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
    add_executable(RANSFieldAudit ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSFieldAudit.cpp)
    target_link_libraries(RANSFieldAudit RANSSpatialDomains ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
    add_executable(RANSWallFieldAudit ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSWallFieldAudit.cpp)
    target_include_directories(RANSWallFieldAudit PRIVATE ${RANS_SOURCE_ROOT}/solvers)
    target_link_libraries(RANSWallFieldAudit RANSKernels RANSSpatialDomains ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
    add_executable(RANSObservationAudit
        ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSObservationAudit.cpp
        ${RANS_SOLVER_ROOT}/EquationSystems/RANSTerrainWall.cpp)
    target_include_directories(RANSObservationAudit PRIVATE ${RANS_SOURCE_ROOT}/solvers)
    target_link_libraries(RANSObservationAudit RANSKernels RANSSpatialDomains
        ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(RANSFieldAudit PRIVATE
            -Wall -Wextra -Wpedantic -Wno-sign-compare -Werror)
        target_compile_options(RANSWallFieldAudit PRIVATE
            -Wall -Wextra -Wpedantic -Wno-sign-compare -Werror)
        target_compile_options(RANSObservationAudit PRIVATE
            -Wall -Wextra -Wpedantic -Wno-sign-compare -Werror)
        set_source_files_properties(
            ${RANS_SOLVER_ROOT}/./EquationSystems/RANSVelocityCorrectionScheme.cpp
            ${RANS_SOLVER_ROOT}/./EquationSystems/RANSAnalyticalWallStrip.cpp
            ${RANS_SOLVER_ROOT}/./EquationSystems/RANSTerrainWall.cpp
            ${RANS_SOLVER_ROOT}/./EquationSystems/RANSMatrixFree.cpp
            ${RANS_SOLVER_ROOT}/./EquationSystems/RANSHistory.cpp
            ${RANS_SOLVER_ROOT}/./EquationSystems/RANSWallFieldLayout.cpp
            ${RANS_SOLVER_ROOT}/./EquationSystems/RANSWallPatch.cpp
            ${RANS_SOLVER_ROOT}/./EquationSystems/RANSWallInterface.cpp
            ${RANS_SOLVER_ROOT}/./EquationSystems/RANSEquilibriumWallSystem.cpp
            ${RANS_SOLVER_ROOT}/./EquationSystems/RANSWallGlobalField.cpp
            ${RANS_SOLVER_ROOT}/./EquationSystems/RANSWallEvolutionPatch.cpp
            ${RANS_SOLVER_ROOT}/./EquationSystems/RANSWallEvolutionEvaluation.cpp
            ${RANS_SOLVER_ROOT}/./EquationSystems/RANSWallEvolutionDiagnostics.cpp
            ${RANS_SOLVER_ROOT}/./EquationSystems/RANSWallEvolutionTransition.cpp
            ${RANS_SOLVER_ROOT}/./EquationSystems/RANSWallEvolutionRobust.cpp
            ${RANS_SOLVER_ROOT}/./EquationSystems/RANSWallEvolutionActions.cpp
            ${RANS_SOLVER_ROOT}/./EquationSystems/RANSWallEvolutionFeasible.cpp
            ${RANS_SOLVER_ROOT}/./EquationSystems/RANSWallEvolutionPicard.cpp
            ${RANS_SOLVER_ROOT}/./EquationSystems/RANSWallEvolutionPimple.cpp
            ${RANS_SOLVER_ROOT}/./EquationSystems/RANSWallEvolutionKBound.cpp
            ${RANS_SOLVER_ROOT}/./EquationSystems/RANSWallEvolutionContact.cpp
            ${RANS_SOLVER_ROOT}/./EquationSystems/RANSWallEvolutionSource.cpp
            ${RANS_SOLVER_ROOT}/./EquationSystems/RANSWallEvolutionReduced.cpp
            ${RANS_SOLVER_ROOT}/./EquationSystems/RANSWallTimeIntegration.cpp
            ${RANS_SOLVER_ROOT}/./EquationSystems/RANSEquilibriumWallEvolution.cpp
            PROPERTIES COMPILE_OPTIONS "-Wall;-Wextra;-Wpedantic;-Wno-sign-compare;-Werror")
    endif()
endif()
enable_testing()
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-evolution-native)
configure_file(${CMAKE_CURRENT_LIST_DIR}/EquilibriumWallFlat.xml
    ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-evolution-native/case.xml COPYONLY)
add_test(NAME RANSEquilibriumTimeNative COMMAND RANSEquilibriumTimeProbe
    --output ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-time.json
    --checkpoint-root ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-time-history)
add_test(NAME RANSEquilibriumEvolutionNative COMMAND RANSEquilibriumEvolutionProbe
    ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-evolution-native/case.xml
    ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-evolution-native/probe.json)
set_tests_properties(RANSEquilibriumTimeNative RANSEquilibriumEvolutionNative PROPERTIES
    PROCESSORS 1 TIMEOUT 300
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-evolution-native
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-native)
configure_file(${CMAKE_CURRENT_LIST_DIR}/EquilibriumWallFlat.xml
    ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-native/case.xml COPYONLY)
add_test(NAME RANSEquilibriumPatchNative COMMAND RANSEquilibriumPatchProbe
    ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-native/case.xml
    ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-native/patch.json)
set_tests_properties(RANSEquilibriumPatchNative PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-native
    PROCESSORS 1 TIMEOUT 300
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
if(RANS_BUILD_SOLVER)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-evolution)
    configure_file(${CMAKE_CURRENT_LIST_DIR}/EquilibriumWallEvolution.xml
        ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-evolution/case.xml COPYONLY)
    add_test(NAME RANSEquilibriumWallEvolutionSmoke COMMAND IncNavierStokesSolver
        ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-evolution/case.xml)
    set_tests_properties(RANSEquilibriumWallEvolutionSmoke PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-evolution
        PROCESSORS 1 TIMEOUT 300
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-stage)
    configure_file(${CMAKE_CURRENT_LIST_DIR}/EquilibriumWallFlat.xml
        ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-stage/case.xml COPYONLY)
    add_test(NAME RANSEquilibriumWallStageSmoke COMMAND IncNavierStokesSolver
        ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-stage/case.xml)
    set_tests_properties(RANSEquilibriumWallStageSmoke PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-stage
        PROCESSORS 1 TIMEOUT 300
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rans-smoke)
    foreach(rans_fixture StandardRANS_MMS CLS_MMS AWFOriginal AWFEquilibrium StandardAWFEquilibrium ABLPreservation ABLPreservation3D StandardResidualMap CLSResidualMap)
        configure_file("${CMAKE_CURRENT_LIST_DIR}/${rans_fixture}.xml"
            "${CMAKE_CURRENT_BINARY_DIR}/rans-smoke/${rans_fixture}.xml" COPYONLY)
    endforeach()
    add_test(NAME RANSStandardTransportSmoke COMMAND IncNavierStokesSolver
        ${CMAKE_CURRENT_BINARY_DIR}/rans-smoke/StandardRANS_MMS.xml)
    set_tests_properties(RANSStandardTransportSmoke PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rans-smoke
        PROCESSORS 1 TIMEOUT 120 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    add_test(NAME RANSCLSTransportSmoke COMMAND IncNavierStokesSolver
        ${CMAKE_CURRENT_BINARY_DIR}/rans-smoke/CLS_MMS.xml)
    add_test(NAME RANSAWFOriginalSmoke COMMAND IncNavierStokesSolver
        ${CMAKE_CURRENT_BINARY_DIR}/rans-smoke/AWFOriginal.xml)
    add_test(NAME RANSAWFEquilibriumSmoke COMMAND IncNavierStokesSolver
        ${CMAKE_CURRENT_BINARY_DIR}/rans-smoke/AWFEquilibrium.xml)
    add_test(NAME RANSStandardAWFEquilibriumSmoke COMMAND IncNavierStokesSolver
        ${CMAKE_CURRENT_BINARY_DIR}/rans-smoke/StandardAWFEquilibrium.xml)
    set_tests_properties(RANSCLSTransportSmoke RANSAWFOriginalSmoke RANSAWFEquilibriumSmoke RANSStandardAWFEquilibriumSmoke PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rans-smoke
        PROCESSORS 1 TIMEOUT 120 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    add_test(NAME RANSABLPreservationSmoke COMMAND IncNavierStokesSolver
        ${CMAKE_CURRENT_BINARY_DIR}/rans-smoke/ABLPreservation.xml)
    add_test(NAME RANSABL3DPerturbationSmoke COMMAND IncNavierStokesSolver
        ${CMAKE_CURRENT_BINARY_DIR}/rans-smoke/ABLPreservation3D.xml)
    set_tests_properties(RANSABLPreservationSmoke RANSABL3DPerturbationSmoke PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rans-smoke
        PROCESSORS 1 TIMEOUT 300 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    add_test(NAME RANSStandardResidualMapSmoke COMMAND IncNavierStokesSolver
        ${CMAKE_CURRENT_BINARY_DIR}/rans-smoke/StandardResidualMap.xml)
    add_test(NAME RANSCLSResidualMapSmoke COMMAND IncNavierStokesSolver
        ${CMAKE_CURRENT_BINARY_DIR}/rans-smoke/CLSResidualMap.xml)
    set_tests_properties(RANSStandardResidualMapSmoke RANSCLSResidualMapSmoke PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rans-smoke
        PROCESSORS 1 TIMEOUT 180 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()
add_test(NAME RANSWallLayerPrototype COMMAND RANSWallLayerProbe
    --output ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-layer.json)
set_tests_properties(RANSWallLayerPrototype PROPERTIES
    PROCESSORS 1 TIMEOUT 60 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
add_test(NAME RANSKernelUnitTests COMMAND RANSKernelProbe
    --output ${CMAKE_CURRENT_BINARY_DIR}/rans-kernels.json)
set_tests_properties(RANSKernelUnitTests PROPERTIES
    PROCESSORS 1 TIMEOUT 120 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
add_test(NAME RANSWallGeometryUnitTests COMMAND RANSWallGeometryProbe
    ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-geometry.json)
set_tests_properties(RANSWallGeometryUnitTests PROPERTIES
    PROCESSORS 1 TIMEOUT 120 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
add_test(NAME RANSObservationKernels COMMAND RANSObservationKernelProbe
    ${CMAKE_CURRENT_BINARY_DIR}/rans-observation-kernels.json)
set_tests_properties(RANSObservationKernels PROPERTIES
    PROCESSORS 1 TIMEOUT 120 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
find_package(Python3 QUIET COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
    add_test(NAME RANSWallTimingChecks COMMAND ${Python3_EXECUTABLE} -B
        ${CMAKE_CURRENT_LIST_DIR}/test_rans_wall_timing.py $<TARGET_FILE:RANSWallTimingChecks>)
    set_tests_properties(RANSWallTimingChecks PROPERTIES PROCESSORS 1 TIMEOUT 30)
endif()
if(Python3_Interpreter_FOUND)
    add_test(NAME RANSAcceptedAWFCases COMMAND ${Python3_EXECUTABLE} -B
        ${CMAKE_CURRENT_LIST_DIR}/test_cls_accepted_awf_cases.py)
    set_tests_properties(RANSAcceptedAWFCases PROPERTIES PROCESSORS 1 TIMEOUT 60
        ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
endif()
add_test(NAME RANSEquilibriumProfilesNative COMMAND RANSEquilibriumProfilesProbe
    --output ${CMAKE_CURRENT_BINARY_DIR}/rans-equilibrium-profiles.json)
set_tests_properties(RANSEquilibriumProfilesNative PROPERTIES
    PROCESSORS 1 TIMEOUT 180 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
if(Python3_Interpreter_FOUND)
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        execute_process(COMMAND ${Python3_EXECUTABLE} -c "import numpy, scipy"
            RESULT_VARIABLE rans_m1_dependencies OUTPUT_QUIET ERROR_QUIET)
        if(rans_m1_dependencies EQUAL 0)
            add_test(NAME RANSCLSEPAWFM1Oracle COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_LIST_DIR}/validate_cls_epawf_m1.py
                --kernel-probe $<TARGET_FILE:RANSKernelProbe>
                --constitutive-probe $<TARGET_FILE:RANSConstitutiveProbe>)
            set_tests_properties(RANSCLSEPAWFM1Oracle PROPERTIES PROCESSORS 1
                TIMEOUT 120 ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
            add_test(NAME RANSCLSEPAWFM2Oracle COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_LIST_DIR}/validate_cls_epawf_m2.py
                --probe $<TARGET_FILE:RANSEquilibriumProfilesProbe>)
            set_tests_properties(RANSCLSEPAWFM2Oracle PROPERTIES PROCESSORS 1
                TIMEOUT 240 ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
        else()
            message(STATUS "CLS EP-AWF M1/M2 oracles require NumPy and SciPy")
        endif()
        if(RANS_BUILD_SOLVER)
            execute_process(COMMAND ${Python3_EXECUTABLE} -c "import numpy, scipy, sympy"
                RESULT_VARIABLE rans_m3_dependencies OUTPUT_QUIET ERROR_QUIET)
            if(rans_m3_dependencies EQUAL 0)
                add_test(NAME RANSCLSEPAWFM8KernelOracle COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_LIST_DIR}/validate_cls_epawf_m8_kernel.py
                    --probe $<TARGET_FILE:RANSRobustVariantProbe>
                    --output ${CMAKE_CURRENT_BINARY_DIR}/cls-epawf-m8-kernel)
                set_tests_properties(RANSCLSEPAWFM8KernelOracle PROPERTIES PROCESSORS 1 TIMEOUT 3600
                    ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
                add_test(NAME RANSCLSEPAWFM8TransverseOracle COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_LIST_DIR}/validate_cls_epawf_m8_transverse.py
                    --build $<TARGET_FILE_DIR:RANSTransverseEvolutionProbe>
                    --output ${CMAKE_CURRENT_BINARY_DIR}/cls-epawf-m8-transverse)
                set_tests_properties(RANSCLSEPAWFM8TransverseOracle PROPERTIES PROCESSORS 1 TIMEOUT 7200
                    ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
                add_test(NAME RANSCLSEPAWFM8NativeOracle COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_LIST_DIR}/validate_cls_epawf_m8.py
                    --build $<TARGET_FILE_DIR:IncNavierStokesSolver> --mode native
                    --output ${CMAKE_CURRENT_BINARY_DIR}/cls-epawf-m8-native)
                set_tests_properties(RANSCLSEPAWFM8NativeOracle PROPERTIES PROCESSORS 2 TIMEOUT 21600
                    ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
                add_test(NAME RANSCLSEPAWFM8RegularizedMMSOracle COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_LIST_DIR}/validate_cls_epawf_m8.py
                    --build $<TARGET_FILE_DIR:IncNavierStokesSolver> --mode mms
                    --variant Regularized --gamma-scale 1
                    --output ${CMAKE_CURRENT_BINARY_DIR}/cls-epawf-m8-regularized-mms)
                set_tests_properties(RANSCLSEPAWFM8RegularizedMMSOracle PROPERTIES PROCESSORS 2 TIMEOUT 21600
                    ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
                add_test(NAME RANSCLSEPAWFM8CappedNativeOracle COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_LIST_DIR}/cls_epawf_m8_capped_native.py
                    --build $<TARGET_FILE_DIR:IncNavierStokesSolver>
                    --output ${CMAKE_CURRENT_BINARY_DIR}/cls-epawf-m8-capped-native)
                set_tests_properties(RANSCLSEPAWFM8CappedNativeOracle PROPERTIES PROCESSORS 1 TIMEOUT 5400
                    ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
                add_test(NAME RANSCLSEPAWFM7NativeOracle COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_LIST_DIR}/validate_cls_epawf_m7.py
                    --build $<TARGET_FILE_DIR:IncNavierStokesSolver>)
                set_tests_properties(RANSCLSEPAWFM7NativeOracle PROPERTIES PROCESSORS 2 TIMEOUT 14400
                    ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
                add_test(NAME RANSCLSEPAWFM7TransverseOracle COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_LIST_DIR}/validate_cls_epawf_m7_transverse.py
                    --build $<TARGET_FILE_DIR:RANSTransverseEvolutionProbe>)
                set_tests_properties(RANSCLSEPAWFM7TransverseOracle PROPERTIES PROCESSORS 1 TIMEOUT 1800
                    ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
                add_test(NAME RANSCLSEPAWFM6KernelOracle COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_LIST_DIR}/validate_cls_epawf_m6_kernel.py
                    --probe $<TARGET_FILE:RANSEquilibriumTransitionProbe>
                    --output ${CMAKE_CURRENT_BINARY_DIR}/cls-epawf-m6-kernel)
                set_tests_properties(RANSCLSEPAWFM6KernelOracle PROPERTIES PROCESSORS 1 TIMEOUT 1800
                    ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
                add_test(NAME RANSCLSEPAWFM6Oracle COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_LIST_DIR}/validate_cls_epawf_m6.py
                    --build $<TARGET_FILE_DIR:RANSCurvedWallProbe>
                    --output ${CMAKE_CURRENT_BINARY_DIR}/cls-epawf-m6)
                set_tests_properties(RANSCLSEPAWFM6Oracle PROPERTIES PROCESSORS 4 TIMEOUT 14400
                    ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
                add_test(NAME RANSCLSEPAWFM6OrientationOracle COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_LIST_DIR}/validate_cls_epawf_m6_orientation.py
                    --build $<TARGET_FILE_DIR:RANSCurvedWallProbe>
                    --output ${CMAKE_CURRENT_BINARY_DIR}/cls-epawf-m6-orientation)
                set_tests_properties(RANSCLSEPAWFM6OrientationOracle PROPERTIES PROCESSORS 1 TIMEOUT 1800
                    ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
                add_test(NAME RANSCLSEPAWFM5NoteOracle COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_LIST_DIR}/validate_cls_epawf_m5_note.py
                    --probe $<TARGET_FILE:RANSCurvedNoteProbe>
                    --profiles ${CMAKE_CURRENT_BINARY_DIR}/rans-equilibrium-profiles.json
                    --output ${CMAKE_CURRENT_BINARY_DIR}/cls-epawf-m5-note)
                set_tests_properties(RANSEquilibriumProfilesNative PROPERTIES FIXTURES_SETUP RANSProfilesForM5)
                set_tests_properties(RANSCLSEPAWFM5NoteOracle PROPERTIES PROCESSORS 1 TIMEOUT 180
                    FIXTURES_REQUIRED RANSProfilesForM5
                    ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
                add_test(NAME RANSCLSEPAWFM5Oracle COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_LIST_DIR}/validate_cls_epawf_m5.py
                    --build $<TARGET_FILE_DIR:RANSCurvedWallProbe>)
                set_tests_properties(RANSCLSEPAWFM5Oracle PROPERTIES PROCESSORS 4
                    TIMEOUT 10800 ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
                add_test(NAME RANSCLSEPAWFM4Oracle COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_LIST_DIR}/validate_cls_epawf_m4.py
                    --build $<TARGET_FILE_DIR:RANSEquilibriumEvolutionProbe>)
                set_tests_properties(RANSCLSEPAWFM4Oracle PROPERTIES PROCESSORS 4
                    TIMEOUT 7200 ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
                add_test(NAME RANSCLSEPAWFM3Oracle COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_LIST_DIR}/validate_cls_epawf_m3.py
                    --build $<TARGET_FILE_DIR:RANSEquilibriumPatchProbe>)
                set_tests_properties(RANSCLSEPAWFM3Oracle PROPERTIES PROCESSORS 1
                    TIMEOUT 1800 ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
            else()
                message(STATUS "CLS EP-AWF M3 oracle requires NumPy, SciPy and SymPy")
            endif()
        endif()
    endif()
    add_test(NAME RANSCaoTamuraCalibrationGuards COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_LIST_DIR}/validate_cao_tamura_calibration.py)
    set_tests_properties(RANSCaoTamuraCalibrationGuards PROPERTIES
        PROCESSORS 1 TIMEOUT 60 ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    add_test(NAME RANSCaoTamuraRouteGuards COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_LIST_DIR}/validate_cao_tamura_routes.py)
    set_tests_properties(RANSCaoTamuraRouteGuards PROPERTIES
        PROCESSORS 1 TIMEOUT 60 ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    add_test(NAME RANSStressTangentOracle COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_LIST_DIR}/validate_rans_stress_tangent.py
        --audit $<TARGET_FILE:RANSStressTangentAudit>)
    set_tests_properties(RANSStressTangentOracle PROPERTIES
        PROCESSORS 1 TIMEOUT 120 ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()
if(RANS_BUILD_SOLVER)
    get_property(rans_all_tests DIRECTORY PROPERTY TESTS)
    set_property(TEST ${rans_all_tests} APPEND PROPERTY ENVIRONMENT
        "LD_LIBRARY_PATH=${RANS_SPATIAL_RUNTIME_DIR}:$ENV{LD_LIBRARY_PATH}")
endif()

# WJ-1 common-kernel probes remain independent of the WJ-2 native pilot.
add_executable(WJKernelProbe ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/WJKernelProbe.cpp)
target_link_libraries(WJKernelProbe RANSKernels)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(WJKernelProbe PRIVATE -Wall -Wextra -Wpedantic -Werror)
    set_source_files_properties(${RANS_SOURCE_ROOT}/library/SolverUtils/Turbulence/WJEARSM.cpp
        PROPERTIES COMPILE_OPTIONS "-Wall;-Wextra;-Wpedantic;-Werror")
endif()
add_test(NAME WJKernelContracts COMMAND WJKernelProbe)
add_executable(WJTerrainGeometryProbe ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/WJTerrainGeometryProbe.cpp)
target_link_libraries(WJTerrainGeometryProbe RANSKernels)
if(RANS_BUILD_SOLVER)
    # Same production consumers as the existing native resolution probe, with
    # only its command-line driver replaced. No installed kernel ABI is mixed.
    get_target_property(wj_h1_native_sources RANSWallResolutionProbe SOURCES)
    list(FILTER wj_h1_native_sources EXCLUDE REGEX "/RANSWallResolutionProbe\\.cpp$")
    add_executable(WJTerrainCoupledProbe
        ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/WJTerrainCoupledProbe.cpp ${wj_h1_native_sources})
    target_include_directories(WJTerrainCoupledProbe PRIVATE ${RANS_SOURCE_ROOT}/solvers)
    target_link_libraries(WJTerrainCoupledProbe RANSKernels RANSSpatialDomains ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
endif()
if(Python3_Interpreter_FOUND)
    execute_process(COMMAND ${Python3_EXECUTABLE} -c "import numpy, scipy, sympy"
        RESULT_VARIABLE wj_dependencies OUTPUT_QUIET ERROR_QUIET)
    if(wj_dependencies EQUAL 0)
        add_test(NAME WJIndependentOracle COMMAND ${Python3_EXECUTABLE} -B
            ${CMAKE_CURRENT_LIST_DIR}/validate_wj1.py --probe $<TARGET_FILE:WJKernelProbe>
            --output ${CMAKE_CURRENT_BINARY_DIR}/wj1-oracle.json)
        set_tests_properties(WJIndependentOracle PROPERTIES PROCESSORS 1 TIMEOUT 240
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
    else()
        message(WARNING "WJIndependentOracle missing required NumPy, SciPy or SymPy")
    endif()
endif()
set_tests_properties(WJKernelContracts PROPERTIES PROCESSORS 1 TIMEOUT 120
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
add_executable(RANSReturnedFieldsProbe ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSReturnedFieldsProbe.cpp)
target_link_libraries(RANSReturnedFieldsProbe RANSKernels)
add_test(NAME RANSReturnedFieldsChecks COMMAND RANSReturnedFieldsProbe)
set_tests_properties(RANSReturnedFieldsChecks PROPERTIES PROCESSORS 1 TIMEOUT 60
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
if(RANS_BUILD_SOLVER AND Python3_Interpreter_FOUND AND wj_dependencies EQUAL 0)
    foreach(wj_policy EP Budget)
        add_test(NAME WJNative${wj_policy} COMMAND ${Python3_EXECUTABLE} -B
            ${CMAKE_CURRENT_LIST_DIR}/test_wj2_native.py
            --solver $<TARGET_FILE:IncNavierStokesSolver> --probe $<TARGET_FILE:RANSWallResolutionProbe>
            --policy ${wj_policy} --output ${CMAKE_CURRENT_BINARY_DIR}/wj2-native-${wj_policy})
        set_tests_properties(WJNative${wj_policy} PROPERTIES PROCESSORS 1 TIMEOUT 1800
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;LD_LIBRARY_PATH=${RANS_SPATIAL_RUNTIME_DIR}:$ENV{LD_LIBRARY_PATH}")
    endforeach()
    add_test(NAME RANSWJSharedLegacyTransaction COMMAND ${Python3_EXECUTABLE} -B
        ${CMAKE_CURRENT_LIST_DIR}/test_wj2_legacy_transaction.py
        --solver $<TARGET_FILE:IncNavierStokesSolver> --probe $<TARGET_FILE:RANSWallResolutionProbe>
        --output ${CMAKE_CURRENT_BINARY_DIR}/wj2-legacy-transaction)
    set_tests_properties(RANSWJSharedLegacyTransaction PROPERTIES PROCESSORS 1 TIMEOUT 1800
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;LD_LIBRARY_PATH=${RANS_SPATIAL_RUNTIME_DIR}:$ENV{LD_LIBRARY_PATH}")
endif()

if(RANS_BUILD_SOLVER AND Python3_Interpreter_FOUND AND wj_dependencies EQUAL 0)
    foreach(wj_policy EP Budget)
        add_test(NAME WJReadOnlyStage${wj_policy} COMMAND ${Python3_EXECUTABLE} -B
            ${CMAKE_CURRENT_LIST_DIR}/run_wj3_native.py
            --solver $<TARGET_FILE:IncNavierStokesSolver> --probe $<TARGET_FILE:RANSWallResolutionProbe>
            --policy ${wj_policy} --mode smoke --output ${CMAKE_CURRENT_BINARY_DIR}/wj3-smoke)
        set_tests_properties(WJReadOnlyStage${wj_policy} PROPERTIES PROCESSORS 1 TIMEOUT 600
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;LD_LIBRARY_PATH=${RANS_SPATIAL_RUNTIME_DIR}:$ENV{LD_LIBRARY_PATH}")
    endforeach()
endif()

if(Python3_Interpreter_FOUND)
    add_test(NAME WJQualificationMetadata COMMAND ${Python3_EXECUTABLE} -B ${CMAKE_CURRENT_LIST_DIR}/test_wj4_metadata.py)
    set_tests_properties(WJQualificationMetadata PROPERTIES PROCESSORS 1 TIMEOUT 30)
endif()
if(Python3_Interpreter_FOUND AND wj_dependencies EQUAL 0)
    add_test(NAME WJLocalNormalMetadata COMMAND ${Python3_EXECUTABLE} -B ${CMAKE_CURRENT_LIST_DIR}/test_wj5_metadata.py)
    set_tests_properties(WJLocalNormalMetadata PROPERTIES PROCESSORS 1 TIMEOUT 60
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
endif()
if(RANS_BUILD_SOLVER AND Python3_Interpreter_FOUND AND wj_dependencies EQUAL 0)
    foreach(wj_policy EP Budget)
        add_test(NAME WJCompositeObservation${wj_policy} COMMAND ${Python3_EXECUTABLE} -B
            ${CMAKE_CURRENT_LIST_DIR}/test_wj4_observations.py
            --solver $<TARGET_FILE:IncNavierStokesSolver> --probe $<TARGET_FILE:RANSWallResolutionProbe>
            --policy ${wj_policy} --output ${CMAKE_CURRENT_BINARY_DIR}/wj4-observations)
        set_tests_properties(WJCompositeObservation${wj_policy} PROPERTIES PROCESSORS 1 TIMEOUT 600
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;LD_LIBRARY_PATH=${RANS_SPATIAL_RUNTIME_DIR}:$ENV{LD_LIBRARY_PATH}")
    endforeach()
    foreach(wj_policy EP Budget)
        add_test(NAME WJLocalNormalObservation${wj_policy} COMMAND ${Python3_EXECUTABLE} -B
            ${CMAKE_CURRENT_LIST_DIR}/test_wj5_observations.py
            --solver $<TARGET_FILE:IncNavierStokesSolver> --probe $<TARGET_FILE:RANSWallResolutionProbe>
            --policy ${wj_policy} --output ${CMAKE_CURRENT_BINARY_DIR}/wj5-observations)
        set_tests_properties(WJLocalNormalObservation${wj_policy} PROPERTIES PROCESSORS 1 TIMEOUT 900
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;LD_LIBRARY_PATH=${RANS_SPATIAL_RUNTIME_DIR}:$ENV{LD_LIBRARY_PATH}")
    endforeach()
    add_test(NAME WJDenseCapacityMetadata COMMAND ${Python3_EXECUTABLE} -B
        ${CMAKE_CURRENT_LIST_DIR}/test_wj6_metadata.py)
    set_tests_properties(WJDenseCapacityMetadata PROPERTIES PROCESSORS 1 TIMEOUT 60
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
    add_test(NAME WJStageLineageContinuationMetadata COMMAND ${Python3_EXECUTABLE} -B
        ${CMAKE_CURRENT_LIST_DIR}/test_wj7_metadata.py)
    set_tests_properties(WJStageLineageContinuationMetadata PROPERTIES PROCESSORS 1 TIMEOUT 90
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
    add_test(NAME WJH1CompositeGeometry COMMAND ${Python3_EXECUTABLE} -B ${CMAKE_CURRENT_LIST_DIR}/validate_wj_h1.py
        --mode geometry --geometry-probe $<TARGET_FILE:WJTerrainGeometryProbe> --output ${CMAKE_CURRENT_BINARY_DIR}/wj-h1-geometry)
    add_test(NAME WJH1RouteGuards COMMAND ${Python3_EXECUTABLE} -B ${CMAKE_CURRENT_LIST_DIR}/validate_wj_h1.py
        --mode guards --native-probe $<TARGET_FILE:WJTerrainCoupledProbe> --output ${CMAKE_CURRENT_BINARY_DIR}/wj-h1-guards)
    add_test(NAME WJH1RuntimeClosure COMMAND ${Python3_EXECUTABLE} -B ${CMAKE_CURRENT_LIST_DIR}/test_wj_h1_runtime.py
        --probe $<TARGET_FILE:WJTerrainCoupledProbe> --library $<TARGET_FILE:RANSSpatialDomains>)
    foreach(wj_h1_kind hill zero plane)
        add_test(NAME WJH1Native_${wj_h1_kind} COMMAND ${Python3_EXECUTABLE} -B ${CMAKE_CURRENT_LIST_DIR}/validate_wj_h1.py
            --mode native --kind ${wj_h1_kind} --native-probe $<TARGET_FILE:WJTerrainCoupledProbe>
            --output ${CMAKE_CURRENT_BINARY_DIR}/wj-h1-${wj_h1_kind})
    endforeach()
    set_tests_properties(WJH1CompositeGeometry WJH1RouteGuards WJH1RuntimeClosure WJH1Native_hill WJH1Native_zero WJH1Native_plane PROPERTIES
        PROCESSORS 1 TIMEOUT 1200 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;LD_LIBRARY_PATH=${RANS_SPATIAL_RUNTIME_DIR}:$ENV{LD_LIBRARY_PATH}")
    foreach(wj_policy EP Budget)
        add_test(NAME WJDenseCapacity${wj_policy} COMMAND ${Python3_EXECUTABLE} -B
            ${CMAKE_CURRENT_LIST_DIR}/test_wj6_capacity.py
            --solver $<TARGET_FILE:IncNavierStokesSolver> --probe $<TARGET_FILE:RANSWallResolutionProbe>
            --policy ${wj_policy} --output ${CMAKE_CURRENT_BINARY_DIR}/wj6-capacity)
        set_tests_properties(WJDenseCapacity${wj_policy} PROPERTIES PROCESSORS 1 TIMEOUT 1800
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;LD_LIBRARY_PATH=${RANS_SPATIAL_RUNTIME_DIR}:$ENV{LD_LIBRARY_PATH}")
    endforeach()
endif()

if(Python3_Interpreter_FOUND)
    add_test(NAME WJH9StartupAudit COMMAND ${Python3_EXECUTABLE} -B ${CMAKE_CURRENT_SOURCE_DIR}/test_wj_h9_startup_audit.py)
    set_tests_properties(WJH9StartupAudit PROPERTIES PROCESSORS 1 TIMEOUT 120
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
    add_test(NAME WJH8LocalTimestep COMMAND ${Python3_EXECUTABLE} -B ${CMAKE_CURRENT_SOURCE_DIR}/test_wj_h8_local_timestep.py)
    set_tests_properties(WJH8LocalTimestep PROPERTIES PROCESSORS 1 TIMEOUT 120
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
    add_test(NAME WJH7NeutralTransport COMMAND ${Python3_EXECUTABLE} -B ${CMAKE_CURRENT_SOURCE_DIR}/test_wj_h7_transport.py)
    set_tests_properties(WJH7NeutralTransport PROPERTIES PROCESSORS 1 TIMEOUT 120
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
    add_test(NAME WJH6NeutralProjection COMMAND ${Python3_EXECUTABLE} -B ${CMAKE_CURRENT_SOURCE_DIR}/test_wj_h6_projection.py)
    set_tests_properties(WJH6NeutralProjection PROPERTIES PROCESSORS 1 TIMEOUT 120
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
    add_test(NAME WJH5NeutralInlet COMMAND ${Python3_EXECUTABLE} -B ${CMAKE_CURRENT_SOURCE_DIR}/test_wj_h5_neutral.py)
    set_tests_properties(WJH5NeutralInlet PROPERTIES PROCESSORS 1 TIMEOUT 120
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
    add_test(NAME WJH4PhysicalOwners COMMAND ${Python3_EXECUTABLE} -B ${CMAKE_CURRENT_SOURCE_DIR}/test_wj_h4_points.py)
    set_tests_properties(WJH4PhysicalOwners PROPERTIES PROCESSORS 1 TIMEOUT 120
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
    add_test(NAME WJH3PressureSeeds COMMAND ${Python3_EXECUTABLE} -B ${CMAKE_CURRENT_SOURCE_DIR}/test_wj_h3_seeds.py)
    set_tests_properties(WJH3PressureSeeds PROPERTIES PROCESSORS 1 TIMEOUT 120
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
    add_test(NAME WJH3PassiveAlgebra COMMAND ${Python3_EXECUTABLE} -B ${CMAKE_CURRENT_SOURCE_DIR}/test_wj_h3_analysis.py)
    set_tests_properties(WJH3PassiveAlgebra PROPERTIES PROCESSORS 1 TIMEOUT 120
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
    add_test(NAME WJH2MeasuredProfiles
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_wj_h2_profiles.py)
    set_tests_properties(WJH2MeasuredProfiles PROPERTIES PROCESSORS 1 TIMEOUT 120
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
endif()
if(RANS_BUILD_SOLVER AND Python3_Interpreter_FOUND)
    add_test(NAME WJH2OpenBoundaries COMMAND ${Python3_EXECUTABLE} -B
        ${CMAKE_CURRENT_SOURCE_DIR}/test_wj_h2_native.py
        --solver $<TARGET_FILE:IncNavierStokesSolver> --probe $<TARGET_FILE:RANSWallResolutionProbe>
        --output ${CMAKE_CURRENT_BINARY_DIR}/wj-h2-open-controls)
    set_tests_properties(WJH2OpenBoundaries PROPERTIES PROCESSORS 1 TIMEOUT 1500
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1;LD_LIBRARY_PATH=${RANS_SPATIAL_RUNTIME_DIR}:$ENV{LD_LIBRARY_PATH}")
endif()
if(Python3_Interpreter_FOUND)
    add_test(NAME WJH2MeasurementContracts COMMAND ${Python3_EXECUTABLE} -B
        ${CMAKE_CURRENT_SOURCE_DIR}/test_wj_h2_measure.py)
    set_tests_properties(WJH2MeasurementContracts PROPERTIES PROCESSORS 1 TIMEOUT 60
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
endif()
if(RANS_BUILD_SOLVER AND Python3_Interpreter_FOUND)
    add_test(NAME WJH2RuntimeClosure COMMAND ${Python3_EXECUTABLE} -B
        ${CMAKE_CURRENT_SOURCE_DIR}/test_wj_h1_runtime.py --probe $<TARGET_FILE:IncNavierStokesSolver>
        --library $<TARGET_FILE:RANSSpatialDomains>)
    set_tests_properties(WJH2RuntimeClosure PROPERTIES PROCESSORS 1 TIMEOUT 60
        ENVIRONMENT "LD_LIBRARY_PATH=${RANS_SPATIAL_RUNTIME_DIR}:$ENV{LD_LIBRARY_PATH}")
endif()

if(Python3_Interpreter_FOUND)
    add_test(NAME WJH10FlatInletInitialization
        COMMAND ${Python3_EXECUTABLE} -B -m unittest -v test_wj_h10_initialization)
    set_tests_properties(WJH10FlatInletInitialization PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROCESSORS 1 TIMEOUT 60
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
endif()

if(Python3_Interpreter_FOUND)
    add_test(NAME WJH11CompatibleContinuation
        COMMAND ${Python3_EXECUTABLE} -B -m unittest -v test_wj_h11_continuation)
    set_tests_properties(WJH11CompatibleContinuation PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROCESSORS 1 TIMEOUT 60
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
endif()

if(Python3_Interpreter_FOUND)
    add_test(NAME WJH12SavedPhaseDiagnosis
        COMMAND ${Python3_EXECUTABLE} -B -m unittest -v test_wj_h12_saved_phase)
    set_tests_properties(WJH12SavedPhaseDiagnosis PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROCESSORS 1 TIMEOUT 60
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
endif()

if(Python3_Interpreter_FOUND)
    add_test(NAME WJH13PairedDomainPreparation
        COMMAND ${Python3_EXECUTABLE} -B -m unittest -v test_wj_h13_domain)
    set_tests_properties(WJH13PairedDomainPreparation PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROCESSORS 1 TIMEOUT 120
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
endif()

if(Python3_Interpreter_FOUND)
    add_test(NAME WJH15PersistentMonitoring
        COMMAND ${Python3_EXECUTABLE} -B -m unittest -v test_wj_h15_monitoring)
    set_tests_properties(WJH15PersistentMonitoring PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROCESSORS 1 TIMEOUT 120
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
    add_test(NAME WJH14FixedHeightExecution
        COMMAND ${Python3_EXECUTABLE} -B -m unittest -v test_wj_h14_execution)
    set_tests_properties(WJH14FixedHeightExecution PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROCESSORS 1 TIMEOUT 120
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;VECLIB_MAXIMUM_THREADS=1")
endif()
