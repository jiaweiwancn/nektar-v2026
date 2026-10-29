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
    ${RANS_SOURCE_ROOT}/library/SolverUtils/Turbulence/AnalyticalWallFunction.cpp)
target_compile_definitions(RANSKernels PRIVATE SOLVER_UTILS_EXPORTS)
add_executable(RANSWallLayerProbe
    ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSWallLayerProbe.cpp)
add_executable(RANSKernelProbe
    ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSKernelProbe.cpp)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    foreach(rans_target RANSKernels RANSWallLayerProbe RANSKernelProbe)
        target_compile_options(${rans_target} PRIVATE
            -Wall -Wextra -Wpedantic -Wno-sign-compare -Werror)
    endforeach()
endif()
target_link_libraries(RANSWallLayerProbe
    RANSKernels ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
target_link_libraries(RANSKernelProbe RANSKernels)
option(RANS_BUILD_SOLVER "Build the complete incompressible solver in isolation" OFF)
if(RANS_BUILD_SOLVER)
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
    target_link_libraries(IncNavierStokesSolver RANSKernels
        ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
    add_executable(RANSFieldAudit ${RANS_SOURCE_ROOT}/library/UnitTests/SolverUtils/RANSFieldAudit.cpp)
    target_link_libraries(RANSFieldAudit ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
    if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
        target_compile_options(RANSFieldAudit PRIVATE
            -Wall -Wextra -Wpedantic -Wno-sign-compare -Werror)
        set_source_files_properties(
            ${RANS_SOLVER_ROOT}/./EquationSystems/RANSVelocityCorrectionScheme.cpp
            PROPERTIES COMPILE_OPTIONS "-Wall;-Wextra;-Wpedantic;-Wno-sign-compare;-Werror")
    endif()
endif()
enable_testing()
if(RANS_BUILD_SOLVER)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rans-smoke)
    add_test(NAME RANSStandardTransportSmoke COMMAND IncNavierStokesSolver
        ${CMAKE_CURRENT_LIST_DIR}/StandardRANS_MMS.xml)
    set_tests_properties(RANSStandardTransportSmoke PROPERTIES
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rans-smoke
        PROCESSORS 1 TIMEOUT 120 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()
add_test(NAME RANSWallLayerPrototype COMMAND RANSWallLayerProbe
    --output ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-layer.json)
set_tests_properties(RANSWallLayerPrototype PROPERTIES
    PROCESSORS 1 TIMEOUT 60 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
add_test(NAME RANSKernelUnitTests COMMAND RANSKernelProbe
    --output ${CMAKE_CURRENT_BINARY_DIR}/rans-kernels.json)
set_tests_properties(RANSKernelUnitTests PROPERTIES
    PROCESSORS 1 TIMEOUT 120 ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
