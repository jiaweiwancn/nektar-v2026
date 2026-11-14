# Small external build against an existing installation. It does not modify
# the shared Nektar++ build or installation while another experiment uses them.
cmake_minimum_required(VERSION 3.12)
project(RANSWallLayerPrototype CXX)
find_package(Nektar++ REQUIRED)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
get_filename_component(RANS_SOURCE_ROOT "${CMAKE_CURRENT_LIST_DIR}/../../../.." ABSOLUTE)
option(RANS_ENABLE_SERIAL_SPARSE "Build the explicitly guarded CT9 serial sparse adapter" OFF)
if(RANS_ENABLE_SERIAL_SPARSE)
    find_package(Eigen3 3.3 REQUIRED NO_MODULE)
    include_directories(SYSTEM ${EIGEN3_INCLUDE_DIRS})
    add_definitions(-DNEKTAR_RANS_HAVE_EIGEN_SPARSE)
endif()
add_executable(RANSReducedModelChecks test_cls_epawf_reduced_model.cpp)
target_include_directories(RANSReducedModelChecks PRIVATE ${RANS_SOURCE_ROOT}/solvers)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(RANSReducedModelChecks PRIVATE -Wall -Wextra -Wpedantic -Wno-sign-compare -Werror)
endif()
add_executable(RANSKBoundQPChecks test_cls_epawf_k_bound_qp.cpp)
target_include_directories(RANSKBoundQPChecks PRIVATE ${RANS_SOURCE_ROOT}/solvers)
add_executable(RANSWallTimingChecks test_rans_wall_timing.cpp)
target_include_directories(RANSWallTimingChecks PRIVATE ${RANS_SOURCE_ROOT}/solvers)
enable_testing()
add_test(NAME RANSKBoundQPChecks COMMAND RANSKBoundQPChecks)
add_test(NAME RANSReducedModelChecks COMMAND RANSReducedModelChecks)
add_definitions(${NEKTAR++_DEFINITIONS} -DENABLE_NEKTAR_EXCEPTIONS)
include_directories(${RANS_SOURCE_ROOT}/library)
include_directories(SYSTEM ${NEKTAR++_INCLUDE_DIRS} ${NEKTAR++_TP_INCLUDE_DIRS})
add_library(RANSKernels STATIC
    ${RANS_SOURCE_ROOT}/library/SolverUtils/Turbulence/RANSModels.cpp
    ${RANS_SOURCE_ROOT}/library/SolverUtils/Turbulence/RANSRobustVariant.cpp
    ${RANS_SOURCE_ROOT}/library/SolverUtils/Turbulence/RANSWallInnerStress.cpp
    ${RANS_SOURCE_ROOT}/library/SolverUtils/Turbulence/AnalyticalWallFunction.cpp
    ${RANS_SOURCE_ROOT}/library/SolverUtils/Turbulence/EquilibriumAWFProfiles.cpp
    ${RANS_SOURCE_ROOT}/library/SolverUtils/Turbulence/EquilibriumAWFTransition.cpp
    ${RANS_SOURCE_ROOT}/library/SolverUtils/Turbulence/RANSWallGeometry.cpp
    ${RANS_SOURCE_ROOT}/library/SolverUtils/Turbulence/RANSWallNormalMap.cpp)
target_compile_definitions(RANSKernels PRIVATE SOLVER_UTILS_EXPORTS)
target_include_directories(RANSKernels PRIVATE ${RANS_SOURCE_ROOT}/solvers)
target_sources(RANSKernels PRIVATE
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallSerialSparse.cpp
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionSparse.cpp)
if(RANS_ENABLE_SERIAL_SPARSE)
    add_executable(CT9SparseAlgebraChecks test_ct_accepted9_sparse.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallFieldLayout.cpp)
    target_include_directories(CT9SparseAlgebraChecks PRIVATE ${RANS_SOURCE_ROOT}/solvers)
    target_link_libraries(CT9SparseAlgebraChecks RANSKernels ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
    add_test(NAME CT9SparseAlgebraChecks COMMAND CT9SparseAlgebraChecks)
    set_tests_properties(CT9SparseAlgebraChecks PROPERTIES PROCESSORS 1 TIMEOUT 60
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;BLIS_NUM_THREADS=1")
endif()
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
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionOpenBoundary.cpp
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
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionOpenBoundary.cpp
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
    ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionOpenBoundary.cpp
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
foreach(rans_r3_probe RANSNonpressureActionProbe RANSWallResolutionProbe RANSPreservationProbe)
    add_executable(${rans_r3_probe}
        ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/${rans_r3_probe}.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallFieldLayout.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallPatch.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallGlobalField.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.cpp
        ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallEvolutionOpenBoundary.cpp
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
target_sources(RANSPreservationProbe PRIVATE
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
        add_library(RANSSpatialDomains SHARED IMPORTED)
        set_target_properties(RANSSpatialDomains PROPERTIES
            IMPORTED_LOCATION "${RANS_PREBUILT_SPATIAL_LIBRARY}")
    else()
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
    # These objects are supplied once by the private RANSKernels archive.
    list(REMOVE_ITEM rans_absolute_sources
        ${RANS_SOLVER_ROOT}/./EquationSystems/RANSWallSerialSparse.cpp
        ${RANS_SOLVER_ROOT}/./EquationSystems/RANSWallEvolutionSparse.cpp)
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
            ${RANS_SOLVER_ROOT}/./EquationSystems/RANSWallEvolutionOpenBoundary.cpp
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
    add_test(NAME RANSImmutableRequestManifest COMMAND ${Python3_EXECUTABLE} -B ${CMAKE_CURRENT_LIST_DIR}/test_rans_request_manifest.py)
    set_tests_properties(RANSImmutableRequestManifest PROPERTIES PROCESSORS 1 TIMEOUT 30)
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
        "LD_LIBRARY_PATH=${CMAKE_CURRENT_BINARY_DIR}/geometry-lib:$ENV{LD_LIBRARY_PATH}")
endif()

# Original A93/A95/E94E95 descriptor, constitutive and wall sampling adapter.
add_executable(RANSNewModelsProbe
    ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSNewModelsProbe.cpp)
target_link_libraries(RANSNewModelsProbe RANSKernels)
add_executable(RANSNativeDiagnosticChecks ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSNativeDiagnosticChecks.cpp)
target_include_directories(RANSNativeDiagnosticChecks PRIVATE ${RANS_SOURCE_ROOT}/solvers)
target_link_libraries(RANSNativeDiagnosticChecks RANSKernels ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
add_test(NAME RANSNativeDiagnosticChecks COMMAND RANSNativeDiagnosticChecks)
set_tests_properties(RANSNativeDiagnosticChecks PROPERTIES PROCESSORS 1 TIMEOUT 30
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(RANSNewModelsProbe PRIVATE -Wall -Wextra -Wpedantic -Wno-sign-compare -Werror)
endif()
find_package(Python3 3.10 REQUIRED COMPONENTS Interpreter)
execute_process(COMMAND ${Python3_EXECUTABLE} -c "import numpy, scipy"
    RESULT_VARIABLE RANS_PYTHON_DEPENDENCY_STATUS ERROR_VARIABLE RANS_PYTHON_DEPENDENCY_ERROR)
if(NOT RANS_PYTHON_DEPENDENCY_STATUS EQUAL 0)
    message(FATAL_ERROR "RANS qualification requires NumPy and SciPy: ${RANS_PYTHON_DEPENDENCY_ERROR}")
endif()
if(Python3_Interpreter_FOUND)
    add_test(NAME RANSPilotEvidenceContract COMMAND ${Python3_EXECUTABLE} -B ${CMAKE_CURRENT_LIST_DIR}/test_rans_pilot_evidence.py)
    add_test(NAME RANSDurationCoverageContract COMMAND ${Python3_EXECUTABLE} -B ${CMAKE_CURRENT_LIST_DIR}/test_rans_duration_contract.py)
    add_test(NAME RANSNonlinearEvidenceContract COMMAND ${Python3_EXECUTABLE} -B ${CMAKE_CURRENT_LIST_DIR}/test_rans_nonlinear_evidence.py)
    add_test(NAME RANSBudgetAccuracyContract COMMAND ${Python3_EXECUTABLE} -B ${CMAKE_CURRENT_LIST_DIR}/test_rans_budget_accuracy_contract.py)
    add_test(NAME RANSOwnershipRequiredScopeContract COMMAND ${Python3_EXECUTABLE} -B ${CMAKE_CURRENT_LIST_DIR}/test_rans_ownership_contract.py)
    add_test(NAME RANSDirectionalRequiredScopeContract COMMAND ${Python3_EXECUTABLE} -B ${CMAKE_CURRENT_LIST_DIR}/test_rans_directional_scope.py)
    set_tests_properties(RANSDirectionalRequiredScopeContract PROPERTIES PROCESSORS 1 TIMEOUT 30)
    set_tests_properties(RANSOwnershipRequiredScopeContract PROPERTIES PROCESSORS 1 TIMEOUT 30
        ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
    set_tests_properties(RANSBudgetAccuracyContract PROPERTIES PROCESSORS 1 TIMEOUT 30
        ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
    set_tests_properties(RANSNonlinearEvidenceContract PROPERTIES PROCESSORS 1 TIMEOUT 30
        ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
    set_tests_properties(RANSDurationCoverageContract PROPERTIES PROCESSORS 1 TIMEOUT 30
        ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
    set_tests_properties(RANSPilotEvidenceContract PROPERTIES PROCESSORS 1 TIMEOUT 30
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
    add_test(NAME RANSQualificationFailureContract COMMAND ${Python3_EXECUTABLE} -B ${CMAKE_CURRENT_LIST_DIR}/test_rans_qualification.py)
    set_tests_properties(RANSQualificationFailureContract PROPERTIES PROCESSORS 1 TIMEOUT 30)
    if(RANS_BUILD_SOLVER)
        add_test(NAME CTAcceptedOpenConservativeNative COMMAND ${Python3_EXECUTABLE} -B
            ${CMAKE_CURRENT_LIST_DIR}/ct_accepted_open.py
            --runtime $<TARGET_FILE_DIR:RANSPreservationProbe>)
        set_tests_properties(CTAcceptedOpenConservativeNative PROPERTIES PROCESSORS 8 TIMEOUT 600
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;PYTHONDONTWRITEBYTECODE=1;TMPDIR=$ENV{TMPDIR};LD_LIBRARY_PATH=$<TARGET_FILE_DIR:RANSSpatialDomains>:$ENV{LD_LIBRARY_PATH}")
        add_test(NAME RANSNativeOwnershipGeometry COMMAND ${Python3_EXECUTABLE} -B
            ${CMAKE_CURRENT_LIST_DIR}/check_rans_native_ownership.py
            --probe $<TARGET_FILE:RANSPreservationProbe> --descriptor-probe $<TARGET_FILE:RANSNewModelsProbe>)
        set_tests_properties(RANSNativeOwnershipGeometry PROPERTIES PROCESSORS 1 TIMEOUT 1800
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;PYTHONDONTWRITEBYTECODE=1;TMPDIR=$ENV{TMPDIR}")
        add_test(NAME RANSDirectionalNativeLayout COMMAND ${Python3_EXECUTABLE} -B
            ${CMAKE_CURRENT_LIST_DIR}/check_rans_directional_layout.py
            --runtime $<TARGET_FILE_DIR:RANSPreservationProbe>)
        set_tests_properties(RANSDirectionalNativeLayout PROPERTIES PROCESSORS 1 TIMEOUT 1200
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;PYTHONDONTWRITEBYTECODE=1;TMPDIR=$ENV{TMPDIR}")
        add_test(NAME RANSForwardInterfaceTolerance COMMAND ${Python3_EXECUTABLE} -B
            ${CMAKE_CURRENT_LIST_DIR}/check_rans_forward_interface_tolerance.py
            --probe $<TARGET_FILE:RANSPreservationProbe> --descriptor-probe $<TARGET_FILE:RANSNewModelsProbe>)
        set_tests_properties(RANSForwardInterfaceTolerance PROPERTIES PROCESSORS 1 TIMEOUT 600
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;PYTHONDONTWRITEBYTECODE=1;TMPDIR=$ENV{TMPDIR}")
        add_test(NAME RANSFrozenBoundaryInterfaceRegression COMMAND ${Python3_EXECUTABLE} -B
            ${CMAKE_CURRENT_LIST_DIR}/test_rans_frozen_boundary_interface.py
            --solver $<TARGET_FILE:IncNavierStokesSolver> --sampler $<TARGET_FILE:RANSPreservationProbe>
            --descriptor-probe $<TARGET_FILE:RANSNewModelsProbe>)
        set_tests_properties(RANSFrozenBoundaryInterfaceRegression PROPERTIES PROCESSORS 1 TIMEOUT 600
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;PYTHONDONTWRITEBYTECODE=1;TMPDIR=$ENV{TMPDIR}")
        add_test(NAME RANSNativeAdmissibilityRegression COMMAND ${Python3_EXECUTABLE} -B
            ${CMAKE_CURRENT_LIST_DIR}/validate_rans_native_audits.py
            --solver $<TARGET_FILE:IncNavierStokesSolver> --descriptor-probe $<TARGET_FILE:RANSNewModelsProbe>)
        set_tests_properties(RANSNativeAdmissibilityRegression PROPERTIES PROCESSORS 1 TIMEOUT 600
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;PYTHONDONTWRITEBYTECODE=1;TMPDIR=$ENV{TMPDIR}")
    endif()
    add_test(NAME RANSNewModelsIndependentOracle COMMAND ${Python3_EXECUTABLE} -B
        ${CMAKE_CURRENT_LIST_DIR}/validate_a93_a95_e94e95.py
        --probe $<TARGET_FILE:RANSNewModelsProbe>
        --output ${CMAKE_CURRENT_BINARY_DIR}/a93-a95-e94e95-kernels.json)
    set_tests_properties(RANSNewModelsIndependentOracle PROPERTIES PROCESSORS 1 TIMEOUT 600
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
endif()

foreach(ct6_native_target RANSEquilibriumEvolutionProbe RANSCurvedWallProbe RANSPressureActionProbe RANSNonpressureActionProbe RANSWallResolutionProbe RANSPreservationProbe)
    if(TARGET ${ct6_native_target})
        target_sources(${ct6_native_target} PRIVATE ${RANS_SOURCE_ROOT}/solvers/IncNavierStokesSolver/EquationSystems/RANSWallManufactured.cpp)
    endif()
endforeach()
add_executable(RANSCurvedAcceptedKernelProbe ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSCurvedAcceptedKernelProbe.cpp)
target_link_libraries(RANSCurvedAcceptedKernelProbe RANSKernels)
add_executable(CT6ManufacturedJetChecks test_ct_accepted6_jet.cpp)
target_include_directories(CT6ManufacturedJetChecks PRIVATE ${RANS_SOURCE_ROOT}/solvers)
add_test(NAME CT6ManufacturedJetChecks COMMAND CT6ManufacturedJetChecks)
if(Python3_Interpreter_FOUND)
    add_test(NAME CT6EvolutionAdmissionChecks COMMAND ${Python3_EXECUTABLE} -B
        ${CMAKE_CURRENT_LIST_DIR}/test_ct_accepted6_gate.py)
    add_test(NAME CT7ScopeAdmissionChecks COMMAND ${Python3_EXECUTABLE} -B
        ${CMAKE_CURRENT_LIST_DIR}/test_ct_accepted7_gate.py)
    add_test(NAME CT8DisplacedAdmissionChecks COMMAND ${Python3_EXECUTABLE} -B
        ${CMAKE_CURRENT_LIST_DIR}/test_ct_accepted8_gate.py)
    set_tests_properties(CT8DisplacedAdmissionChecks PROPERTIES PROCESSORS 1 TIMEOUT 60
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;BLIS_NUM_THREADS=1;PYTHONDONTWRITEBYTECODE=1")
    set_tests_properties(CT7ScopeAdmissionChecks PROPERTIES PROCESSORS 1 TIMEOUT 60
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;NUMEXPR_NUM_THREADS=1;BLIS_NUM_THREADS=1;PYTHONDONTWRITEBYTECODE=1")
    set_tests_properties(CT6EvolutionAdmissionChecks PROPERTIES PROCESSORS 1 TIMEOUT 60
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1;PYTHONDONTWRITEBYTECODE=1")
endif()
