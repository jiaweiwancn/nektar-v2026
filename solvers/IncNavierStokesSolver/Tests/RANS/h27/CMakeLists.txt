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
add_executable(WJH27PrincipalProbe principal.cpp)
target_link_libraries(WJH27PrincipalProbe PRIVATE WJNativeFlow)
target_compile_options(WJH27PrincipalProbe PRIVATE
    -O3 -ffp-contract=off -fno-fast-math -Wall -Wextra)
