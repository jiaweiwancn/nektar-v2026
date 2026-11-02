# Small external build against an existing installation. It does not modify
# the shared Nektar++ build or installation while another experiment uses them.
cmake_minimum_required(VERSION 3.12)
project(RANSWallLayerPrototype CXX)
find_package(Nektar++ REQUIRED)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
get_filename_component(RANS_SOURCE_ROOT "${CMAKE_CURRENT_LIST_DIR}/../../../.." ABSOLUTE)
add_definitions(${NEKTAR++_DEFINITIONS} -DENABLE_NEKTAR_EXCEPTIONS)
include_directories(${RANS_SOURCE_ROOT}/library)
include_directories(SYSTEM ${NEKTAR++_INCLUDE_DIRS} ${NEKTAR++_TP_INCLUDE_DIRS})
add_library(RANSKernels STATIC
    ${RANS_SOURCE_ROOT}/library/SolverUtils/Turbulence/RANSModels.cpp
    ${RANS_SOURCE_ROOT}/library/SolverUtils/Turbulence/AnalyticalWallFunction.cpp
    ${RANS_SOURCE_ROOT}/library/SolverUtils/Turbulence/EquilibriumAWFProfiles.cpp
    ${RANS_SOURCE_ROOT}/library/SolverUtils/Turbulence/RANSWallGeometry.cpp)
target_compile_definitions(RANSKernels PRIVATE SOLVER_UTILS_EXPORTS)
add_executable(RANSWallLayerProbe
    ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSWallLayerProbe.cpp)
add_executable(RANSKernelProbe
    ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSKernelProbe.cpp)
add_executable(RANSWallGeometryProbe
    ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSWallGeometryProbe.cpp)
add_executable(RANSObservationKernelProbe
    ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSObservationKernelProbe.cpp)
add_executable(RANSStressTangentAudit
    ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSStressTangentAudit.cpp)
add_executable(RANSConstitutiveProbe
    ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSConstitutiveProbe.cpp)
add_executable(RANSEquilibriumProfilesProbe
    ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSEquilibriumProfilesProbe.cpp)
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
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    foreach(rans_target RANSKernels RANSWallLayerProbe RANSKernelProbe RANSWallGeometryProbe RANSObservationKernelProbe RANSStressTangentAudit RANSConstitutiveProbe RANSEquilibriumProfilesProbe RANSEquilibriumPatchProbe RANSEquilibriumTimeProbe)
        target_compile_options(${rans_target} PRIVATE
            -Wall -Wextra -Wpedantic -Wno-sign-compare -Werror)
    endforeach()
endif()
target_link_libraries(RANSWallLayerProbe
    RANSKernels ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
target_link_libraries(RANSKernelProbe RANSKernels)
target_link_libraries(RANSWallGeometryProbe RANSKernels)
target_link_libraries(RANSObservationKernelProbe RANSKernels)
target_link_libraries(RANSStressTangentAudit RANSKernels)
target_link_libraries(RANSConstitutiveProbe RANSKernels)
target_link_libraries(RANSEquilibriumProfilesProbe RANSKernels ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
option(RANS_BUILD_SOLVER "Build the complete incompressible solver in isolation" OFF)
if(RANS_BUILD_SOLVER)
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
    # Imported SpatialDomains exposes its partitioner dependencies. Retain those
    # while replacing only that library, rather than loading two copies.
    get_target_property(rans_spatial_dependencies SpatialDomains INTERFACE_LINK_LIBRARIES)
    target_link_libraries(RANSSpatialDomains StdRegions ${rans_spatial_dependencies}
        ${NEKTAR++_TP_LIBRARIES})
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
            ${RANS_SOLVER_ROOT}/./EquationSystems/RANSWallTimeIntegration.cpp
            ${RANS_SOLVER_ROOT}/./EquationSystems/RANSEquilibriumWallEvolution.cpp
            PROPERTIES COMPILE_OPTIONS "-Wall;-Wextra;-Wpedantic;-Wno-sign-compare;-Werror")
    endif()
endif()
enable_testing()
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
