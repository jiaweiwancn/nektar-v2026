IF(TARGET SolverUtils)
    ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS)
    ADD_NEKTAR_EXECUTABLE(SolverUtilsUnitTests
        COMPONENT unit-test DEPENDS SolverUtils
        SOURCES main.cpp TestHIT.cpp)
    TARGET_COMPILE_DEFINITIONS(SolverUtilsUnitTests PRIVATE
        HIT_TEST_FIXTURE_DIR="${CMAKE_SOURCE_DIR}/solvers/IncNavierStokesSolver/Tests/HITInflow")
    ADD_TEST(NAME HITCoreUnitTests
        COMMAND SolverUtilsUnitTests --run_test=HITCore --detect_memory_leaks=0)

    # Exercise the actual boundary objects without requiring solver libraries.
    SET(HIT_BOUNDARY_DIR ${CMAKE_SOURCE_DIR}/solvers/IncNavierStokesSolver/BoundaryConditions)
    ADD_NEKTAR_EXECUTABLE(HITBoundaryUnitTests
        COMPONENT unit-test DEPENDS SolverUtils
        SOURCES main.cpp TestHITBoundary.cpp
        ${HIT_BOUNDARY_DIR}/IncBaseCondition.cpp
        ${HIT_BOUNDARY_DIR}/IncBoundaryConditions.cpp
        ${HIT_BOUNDARY_DIR}/HITInflow.cpp
        ${HIT_BOUNDARY_DIR}/HITInflowConfig.cpp
        ${HIT_BOUNDARY_DIR}/HITInflowGeometry.cpp)
    TARGET_INCLUDE_DIRECTORIES(HITBoundaryUnitTests PRIVATE ${CMAKE_SOURCE_DIR}/solvers)
    TARGET_COMPILE_DEFINITIONS(HITBoundaryUnitTests PRIVATE
        HIT_TEST_FIXTURE_DIR="${CMAKE_SOURCE_DIR}/solvers/IncNavierStokesSolver/Tests/HITInflow")
    ADD_TEST(NAME HITBoundaryUnitTests
        COMMAND HITBoundaryUnitTests --detect_memory_leaks=0)
ENDIF()
