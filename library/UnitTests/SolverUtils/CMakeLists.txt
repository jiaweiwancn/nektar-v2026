IF(TARGET SolverUtils)
    ADD_DEFINITIONS(-DENABLE_NEKTAR_EXCEPTIONS)
    ADD_NEKTAR_EXECUTABLE(SolverUtilsUnitTests
        COMPONENT unit-test DEPENDS SolverUtils
        SOURCES main.cpp TestHIT.cpp)
    TARGET_COMPILE_DEFINITIONS(SolverUtilsUnitTests PRIVATE
        HIT_TEST_FIXTURE_DIR="${CMAKE_SOURCE_DIR}/solvers/IncNavierStokesSolver/Tests/HITInflow")
    ADD_TEST(NAME HITCoreUnitTests
        COMMAND SolverUtilsUnitTests --run_test=HITCore --detect_memory_leaks=0)
ENDIF()
