CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_INCNAVIERSTOKES
    "Build the Incompressible Navier-Stokes solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF (NOT NEKTAR_BUILD_SOLVERS)
    SET(NEKTAR_SOLVER_INCNAVIERSTOKES OFF CACHE INTERNAL "")
ENDIF()

IF( NEKTAR_SOLVER_INCNAVIERSTOKES )
    ADD_SOLVER_EXECUTABLE(IncNavierStokesSolver
        SOURCES ./IncNavierStokesSolver.cpp
        LIBRARY_SOURCES
        ./EquationSystems/CoupledLinearNS.cpp
        ./EquationSystems/CoupledLocalToGlobalC0ContMap.cpp
        ./EquationSystems/IncNavierStokes.cpp
        ./EquationSystems/VelocityCorrectionScheme.cpp
        ./EquationSystems/VelocityCorrectionSchemeImplicit.cpp
        ./EquationSystems/VelocityCorrectionSchemeWeakPressure.cpp
        ./EquationSystems/RANSVelocityCorrectionScheme.cpp
        ./EquationSystems/RANSAnalyticalWallStrip.cpp
        ./EquationSystems/RANSTerrainWall.cpp
        ./EquationSystems/RANSMatrixFree.cpp
        ./EquationSystems/RANSHistory.cpp
        ./EquationSystems/RANSWallFieldLayout.cpp
        ./EquationSystems/RANSWallPatch.cpp
        ./EquationSystems/RANSWallInterface.cpp
        ./EquationSystems/RANSEquilibriumWallSystem.cpp
        ./EquationSystems/RANSWallGlobalField.cpp
        ./EquationSystems/RANSWallEvolutionPatch.cpp
        ./EquationSystems/RANSWallEvolutionDiagnostics.cpp
        ./EquationSystems/RANSWallEvolutionTransition.cpp
        ./EquationSystems/RANSWallEvolutionRobust.cpp
        ./EquationSystems/RANSWallEvolutionActions.cpp
        ./EquationSystems/RANSWallEvolutionFeasible.cpp
        ./EquationSystems/RANSWallTimeIntegration.cpp
        ./EquationSystems/RANSEquilibriumWallEvolution.cpp
        ./EquationSystems/SmoothedProfileMethod.cpp
        ./EquationSystems/VCSMapping.cpp
        ./EquationSystems/Extrapolate.cpp
        ./EquationSystems/StandardExtrapolate.cpp
        ./EquationSystems/MappingExtrapolate.cpp
        ./EquationSystems/SubSteppingExtrapolate.cpp
        ./EquationSystems/SubSteppingExtrapolateWeakPressure.cpp
        ./EquationSystems/WeakPressureExtrapolate.cpp
        ./EquationSystems/ImplicitExtrapolate.cpp
        ./AdvectionTerms/AdjointAdvection.cpp
        ./AdvectionTerms/LinearisedAdvection.cpp
        ./AdvectionTerms/NavierStokesAdvection.cpp
        ./AdvectionTerms/SkewSymmetricAdvection.cpp
        ./AdvectionTerms/AlternateSkewAdvection.cpp
        ./AdvectionTerms/NoAdvection.cpp
        ./Filters/FilterReynoldsStresses.cpp
        ./Filters/FilterMovingBody.cpp
        ./Filters/FilterAeroForcesSPM.cpp
        ./Forcing/ForcingIncNSSyntheticEddy.cpp
        ./Forcing/ForcingMovingBody.cpp
        ./Forcing/ForcingStabilityCoupledLNS.cpp
        ./BoundaryConditions/IncBoundaryConditions.cpp
        ./BoundaryConditions/IncBaseCondition.cpp
        ./BoundaryConditions/HITInflow.cpp
        ./BoundaryConditions/HITInflowConfig.cpp
        ./BoundaryConditions/HITInflowGeometry.cpp
        ./BoundaryConditions/HITInflowRestart.cpp
        ./BoundaryConditions/MovingFrameFar.cpp
        ./BoundaryConditions/StaticWall.cpp
        ./BoundaryConditions/TransMovingWall.cpp
        ./BoundaryConditions/MovingFrameWall.cpp
        )

    IF(NEKTAR_BUILD_SOLVER_LIBS)
        SET(IncNavierStokesSolverHeaders
            AdvectionTerms/AdjointAdvection.h
            AdvectionTerms/AlternateSkewAdvection.h
            AdvectionTerms/LinearisedAdvection.h
            AdvectionTerms/NavierStokesAdvection.h
            AdvectionTerms/NoAdvection.h
            AdvectionTerms/SkewSymmetricAdvection.h
            EquationSystems/CoupledLinearNS.h
            EquationSystems/CoupledLocalToGlobalC0ContMap.h
            EquationSystems/Extrapolate.h
            EquationSystems/ImplicitExtrapolate.h
            EquationSystems/IncNavierStokes.h
            EquationSystems/MappingExtrapolate.h
            EquationSystems/SmoothedProfileMethod.h
            EquationSystems/StandardExtrapolate.h
            EquationSystems/SubSteppingExtrapolate.h
            EquationSystems/SubSteppingExtrapolateWeakPressure.h
            EquationSystems/VCSMapping.h
            EquationSystems/VelocityCorrectionScheme.h
            EquationSystems/VelocityCorrectionSchemeImplicit.h
            EquationSystems/VelocityCorrectionSchemeWeakPressure.h
            EquationSystems/RANSVelocityCorrectionScheme.h
            EquationSystems/RANSAnalyticalWallStrip.h
            EquationSystems/RANSMatrixFree.h
            EquationSystems/RANSHistory.h
            EquationSystems/RANSEquilibriumWallSystem.h
            EquationSystems/RANSWallFieldLayout.h
            EquationSystems/RANSWallPatch.h
            EquationSystems/RANSWallInterface.h
            EquationSystems/RANSWallGlobalField.h
            EquationSystems/RANSWallEvolutionPatch.h
            EquationSystems/RANSWallTimeIntegration.h
            EquationSystems/RANSEquilibriumWallEvolution.h
            EquationSystems/WeakPressureExtrapolate.h
            Filters/FilterAeroForcesSPM.h
            Filters/FilterMovingBody.h
            Filters/FilterReynoldsStresses.h
            Forcing/ForcingIncNSSyntheticEddy.h
            Forcing/ForcingMovingBody.h
            Forcing/ForcingStabilityCoupledLNS.h
            BoundaryConditions/IncBoundaryConditions.h
            BoundaryConditions/IncBaseCondition.h
            BoundaryConditions/HITInflow.h
            BoundaryConditions/HITInflowConfig.h
            BoundaryConditions/HITInflowGeometry.h
            BoundaryConditions/HITInflowRestart.h
            BoundaryConditions/MovingFrameFar.h
            BoundaryConditions/StaticWall.h
            BoundaryConditions/TransMovingWall.h
            BoundaryConditions/MovingFrameWall.h)
        # Create separate library. No need to pass source files, as the executable's
        # LIBRARY_SOURCES objects will be reused
        ADD_SOLVER_LIBRARY(IncNavierStokesSolver HEADERS
                          ${IncNavierStokesSolverHeaders} DEPENDS SolverUtils)
    ENDIF()

    SUBDIRS(Utilities)

    # The following test is disabled pending better integration with the CI system
    ADD_NEKTAR_PERFORMANCE_TEST(Perf_ChanFlow_3DH1D_pRef)
    ADD_NEKTAR_PERFORMANCE_TEST(Perf_ChanStability)

    ADD_NEKTAR_TEST(Kovasznay_Flow_3modes)
    ADD_NEKTAR_TEST(HITInflowMeanOnly3D)
    ADD_NEKTAR_TEST(HITInflowManufactured3D)
    ADD_NEKTAR_TEST(HITInflowHomogeneous1DMean)

    # Native regression tools also build when the optional unit tests are off.
    IF(NEKTAR_BUILD_TESTS)
        ADD_NEKTAR_EXECUTABLE(RANSWallLayerProbe COMPONENT tests DEPENDS SolverUtils
            SOURCES ${CMAKE_SOURCE_DIR}/library/UnitTests/SolverUtils/RANSWallLayerProbe.cpp)
        TARGET_COMPILE_DEFINITIONS(RANSWallLayerProbe PRIVATE ENABLE_NEKTAR_EXCEPTIONS)
        ADD_TEST(NAME RANSWallLayerPrototype COMMAND RANSWallLayerProbe
            --output ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-layer.json)
        SET_TESTS_PROPERTIES(RANSWallLayerPrototype PROPERTIES
            PROCESSORS 1 TIMEOUT 60
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
        ADD_NEKTAR_EXECUTABLE(RANSKernelProbe COMPONENT tests DEPENDS SolverUtils
            SOURCES ${CMAKE_SOURCE_DIR}/library/UnitTests/SolverUtils/RANSKernelProbe.cpp)
        TARGET_COMPILE_DEFINITIONS(RANSKernelProbe PRIVATE ENABLE_NEKTAR_EXCEPTIONS)
        ADD_TEST(NAME RANSKernelUnitTests COMMAND RANSKernelProbe
            --output ${CMAKE_CURRENT_BINARY_DIR}/rans-kernels.json)
        SET_TESTS_PROPERTIES(RANSKernelUnitTests PROPERTIES
            PROCESSORS 1 TIMEOUT 120
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
        ADD_NEKTAR_EXECUTABLE(RANSRobustVariantProbe COMPONENT tests DEPENDS SolverUtils
            SOURCES ${CMAKE_SOURCE_DIR}/library/UnitTests/SolverUtils/RANSRobustVariantProbe.cpp)
        TARGET_COMPILE_DEFINITIONS(RANSRobustVariantProbe PRIVATE ENABLE_NEKTAR_EXCEPTIONS)
        ADD_NEKTAR_EXECUTABLE(RANSConstitutiveProbe COMPONENT tests DEPENDS SolverUtils
            SOURCES ${CMAKE_SOURCE_DIR}/library/UnitTests/SolverUtils/RANSConstitutiveProbe.cpp)
        TARGET_COMPILE_DEFINITIONS(RANSConstitutiveProbe PRIVATE ENABLE_NEKTAR_EXCEPTIONS)
        ADD_NEKTAR_EXECUTABLE(RANSEquilibriumProfilesProbe COMPONENT tests DEPENDS SolverUtils
            SOURCES ${CMAKE_SOURCE_DIR}/library/UnitTests/SolverUtils/RANSEquilibriumProfilesProbe.cpp)
        TARGET_COMPILE_DEFINITIONS(RANSEquilibriumProfilesProbe PRIVATE ENABLE_NEKTAR_EXCEPTIONS)
        ADD_TEST(NAME RANSEquilibriumProfilesNative COMMAND RANSEquilibriumProfilesProbe
            --output ${CMAKE_CURRENT_BINARY_DIR}/rans-equilibrium-profiles.json)
        SET_TESTS_PROPERTIES(RANSEquilibriumProfilesNative PROPERTIES PROCESSORS 1 TIMEOUT 180
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
        ADD_NEKTAR_EXECUTABLE(RANSEquilibriumPatchProbe COMPONENT tests DEPENDS SolverUtils
            SOURCES ${CMAKE_SOURCE_DIR}/library/UnitTests/SolverUtils/RANSEquilibriumPatchProbe.cpp
                EquationSystems/RANSWallFieldLayout.cpp
                EquationSystems/RANSWallPatch.cpp
                EquationSystems/RANSWallInterface.cpp)
        TARGET_INCLUDE_DIRECTORIES(RANSEquilibriumPatchProbe PRIVATE ${CMAKE_SOURCE_DIR}/solvers)
        TARGET_COMPILE_DEFINITIONS(RANSEquilibriumPatchProbe PRIVATE ENABLE_NEKTAR_EXCEPTIONS)
        ADD_NEKTAR_EXECUTABLE(RANSEquilibriumTimeProbe COMPONENT tests DEPENDS SolverUtils
            SOURCES ${CMAKE_SOURCE_DIR}/library/UnitTests/SolverUtils/RANSEquilibriumTimeProbe.cpp
                EquationSystems/RANSWallFieldLayout.cpp
                EquationSystems/RANSWallTimeIntegration.cpp
                EquationSystems/RANSHistory.cpp)
        ADD_NEKTAR_EXECUTABLE(RANSTransverseEvolutionProbe COMPONENT tests DEPENDS SolverUtils
            SOURCES ${CMAKE_SOURCE_DIR}/library/UnitTests/SolverUtils/RANSTransverseEvolutionProbe.cpp
                EquationSystems/RANSWallFieldLayout.cpp
                EquationSystems/RANSWallTimeIntegration.cpp
                EquationSystems/RANSHistory.cpp)
        ADD_NEKTAR_EXECUTABLE(RANSEquilibriumEvolutionProbe COMPONENT tests DEPENDS SolverUtils
            SOURCES ${CMAKE_SOURCE_DIR}/library/UnitTests/SolverUtils/RANSEquilibriumEvolutionProbe.cpp
                EquationSystems/RANSWallFieldLayout.cpp
                EquationSystems/RANSWallPatch.cpp
                EquationSystems/RANSWallGlobalField.cpp
                EquationSystems/RANSWallEvolutionPatch.cpp
                EquationSystems/RANSWallEvolutionTransition.cpp
                EquationSystems/RANSWallEvolutionRobust.cpp)
        TARGET_SOURCES(RANSEquilibriumEvolutionProbe PRIVATE
            EquationSystems/RANSWallEvolutionActions.cpp
            EquationSystems/RANSWallEvolutionFeasible.cpp
            EquationSystems/RANSWallEvolutionDiagnostics.cpp)
        FOREACH(rans_m4_target RANSEquilibriumTimeProbe RANSEquilibriumEvolutionProbe RANSTransverseEvolutionProbe)
            TARGET_INCLUDE_DIRECTORIES(${rans_m4_target} PRIVATE ${CMAKE_SOURCE_DIR}/solvers)
            TARGET_COMPILE_DEFINITIONS(${rans_m4_target} PRIVATE ENABLE_NEKTAR_EXCEPTIONS)
        ENDFOREACH()
        FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-evolution-native)
        CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/Tests/RANS/EquilibriumWallFlat.xml
            ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-evolution-native/case.xml COPYONLY)
        ADD_TEST(NAME RANSEquilibriumTimeNative COMMAND RANSEquilibriumTimeProbe
            --output ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-time.json
            --checkpoint-root ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-time-history)
        ADD_TEST(NAME RANSEquilibriumEvolutionNative COMMAND RANSEquilibriumEvolutionProbe
            ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-evolution-native/case.xml
            ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-evolution-native/probe.json)
        SET_TESTS_PROPERTIES(RANSEquilibriumTimeNative RANSEquilibriumEvolutionNative PROPERTIES
            PROCESSORS 1 TIMEOUT 300
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-evolution-native
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
        FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-evolution)
        CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/Tests/RANS/EquilibriumWallEvolution.xml
            ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-evolution/case.xml COPYONLY)
        ADD_TEST(NAME RANSEquilibriumWallEvolutionSmoke COMMAND IncNavierStokesSolver
            ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-evolution/case.xml)
        SET_TESTS_PROPERTIES(RANSEquilibriumWallEvolutionSmoke PROPERTIES
            PROCESSORS 1 TIMEOUT 300
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-evolution
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
        FOREACH(rans_wall_case native stage)
            FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-${rans_wall_case})
            CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/Tests/RANS/EquilibriumWallFlat.xml
                ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-${rans_wall_case}/case.xml COPYONLY)
        ENDFOREACH()
        ADD_TEST(NAME RANSEquilibriumPatchNative COMMAND RANSEquilibriumPatchProbe
            ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-native/case.xml
            ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-native/patch.json)
        SET_TESTS_PROPERTIES(RANSEquilibriumPatchNative PROPERTIES PROCESSORS 1 TIMEOUT 300
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-native
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
        ADD_TEST(NAME RANSEquilibriumWallStageSmoke COMMAND IncNavierStokesSolver
            ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-stage/case.xml)
        SET_TESTS_PROPERTIES(RANSEquilibriumWallStageSmoke PROPERTIES PROCESSORS 1 TIMEOUT 300
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rans-wall-stage
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
        FIND_PACKAGE(Python3 QUIET COMPONENTS Interpreter)
        IF(Python3_Interpreter_FOUND AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
            EXECUTE_PROCESS(COMMAND ${Python3_EXECUTABLE} -c "import numpy, scipy"
                RESULT_VARIABLE RANS_M1_DEPENDENCIES OUTPUT_QUIET ERROR_QUIET)
            IF(RANS_M1_DEPENDENCIES EQUAL 0)
                ADD_TEST(NAME RANSCLSEPAWFM1Oracle COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/Tests/RANS/validate_cls_epawf_m1.py
                    --kernel-probe $<TARGET_FILE:RANSKernelProbe>
                    --constitutive-probe $<TARGET_FILE:RANSConstitutiveProbe>)
                SET_TESTS_PROPERTIES(RANSCLSEPAWFM1Oracle PROPERTIES PROCESSORS 1
                    TIMEOUT 120 ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
                ADD_TEST(NAME RANSCLSEPAWFM2Oracle COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/Tests/RANS/validate_cls_epawf_m2.py
                    --probe $<TARGET_FILE:RANSEquilibriumProfilesProbe>)
                SET_TESTS_PROPERTIES(RANSCLSEPAWFM2Oracle PROPERTIES PROCESSORS 1
                    TIMEOUT 240 ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
            ELSE()
                MESSAGE(STATUS "CLS EP-AWF M1/M2 oracles require NumPy and SciPy")
            ENDIF()
            EXECUTE_PROCESS(COMMAND ${Python3_EXECUTABLE} -c "import numpy, scipy, sympy"
                RESULT_VARIABLE RANS_M3_DEPENDENCIES OUTPUT_QUIET ERROR_QUIET)
            IF(RANS_M3_DEPENDENCIES EQUAL 0)
                ADD_TEST(NAME RANSCLSEPAWFM4Oracle COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/Tests/RANS/validate_cls_epawf_m4.py
                    --build $<TARGET_FILE_DIR:RANSEquilibriumEvolutionProbe>)
                SET_TESTS_PROPERTIES(RANSCLSEPAWFM4Oracle PROPERTIES PROCESSORS 4
                    TIMEOUT 7200 ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
                ADD_TEST(NAME RANSCLSEPAWFM3Oracle COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/Tests/RANS/validate_cls_epawf_m3.py
                    --build $<TARGET_FILE_DIR:RANSEquilibriumPatchProbe>)
                SET_TESTS_PROPERTIES(RANSCLSEPAWFM3Oracle PROPERTIES PROCESSORS 1
                    TIMEOUT 1800 ENVIRONMENT "PYTHONDONTWRITEBYTECODE=1;OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1;BLIS_NUM_THREADS=1")
            ELSE()
                MESSAGE(STATUS "CLS EP-AWF M3 oracle requires NumPy, SciPy and SymPy")
            ENDIF()
        ENDIF()
        ADD_NEKTAR_EXECUTABLE(RANSFieldAudit COMPONENT tests DEPENDS MultiRegions
            SOURCES ${CMAKE_SOURCE_DIR}/library/UnitTests/SolverUtils/RANSFieldAudit.cpp)
        FILE(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rans-smoke)
        FOREACH(rans_fixture StandardRANS_MMS CLS_MMS AWFOriginal AWFEquilibrium StandardAWFEquilibrium ABLPreservation ABLPreservation3D)
            CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/Tests/RANS/${rans_fixture}.xml"
                "${CMAKE_CURRENT_BINARY_DIR}/rans-smoke/${rans_fixture}.xml" COPYONLY)
        ENDFOREACH()
        ADD_TEST(NAME RANSStandardTransportSmoke COMMAND IncNavierStokesSolver
            ${CMAKE_CURRENT_BINARY_DIR}/rans-smoke/StandardRANS_MMS.xml)
        SET_TESTS_PROPERTIES(RANSStandardTransportSmoke PROPERTIES
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rans-smoke
            PROCESSORS 1 TIMEOUT 120
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
        ADD_TEST(NAME RANSCLSTransportSmoke COMMAND IncNavierStokesSolver
            ${CMAKE_CURRENT_BINARY_DIR}/rans-smoke/CLS_MMS.xml)
        ADD_TEST(NAME RANSAWFOriginalSmoke COMMAND IncNavierStokesSolver
            ${CMAKE_CURRENT_BINARY_DIR}/rans-smoke/AWFOriginal.xml)
        ADD_TEST(NAME RANSAWFEquilibriumSmoke COMMAND IncNavierStokesSolver
            ${CMAKE_CURRENT_BINARY_DIR}/rans-smoke/AWFEquilibrium.xml)
        ADD_TEST(NAME RANSStandardAWFEquilibriumSmoke COMMAND IncNavierStokesSolver
            ${CMAKE_CURRENT_BINARY_DIR}/rans-smoke/StandardAWFEquilibrium.xml)
        SET_TESTS_PROPERTIES(RANSCLSTransportSmoke RANSAWFOriginalSmoke RANSAWFEquilibriumSmoke RANSStandardAWFEquilibriumSmoke PROPERTIES
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rans-smoke
            PROCESSORS 1 TIMEOUT 120
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
        ADD_TEST(NAME RANSABLPreservationSmoke COMMAND IncNavierStokesSolver
            ${CMAKE_CURRENT_BINARY_DIR}/rans-smoke/ABLPreservation.xml)
        ADD_TEST(NAME RANSABL3DPerturbationSmoke COMMAND IncNavierStokesSolver
            ${CMAKE_CURRENT_BINARY_DIR}/rans-smoke/ABLPreservation3D.xml)
        SET_TESTS_PROPERTIES(RANSABLPreservationSmoke RANSABL3DPerturbationSmoke PROPERTIES
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/rans-smoke
            PROCESSORS 1 TIMEOUT 300
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
        SET(HIT_TEST_DIR ${CMAKE_SOURCE_DIR}/library/UnitTests/SolverUtils)
        SET(HIT_TEST_BOUNDARY_SOURCES
            BoundaryConditions/IncBaseCondition.cpp
            BoundaryConditions/IncBoundaryConditions.cpp
            BoundaryConditions/HITInflow.cpp
            BoundaryConditions/HITInflowConfig.cpp
            BoundaryConditions/HITInflowGeometry.cpp
            BoundaryConditions/HITInflowRestart.cpp)
        ADD_NEKTAR_EXECUTABLE(HITBoundaryProbe COMPONENT tests DEPENDS SolverUtils
            SOURCES ${HIT_TEST_DIR}/HITBoundaryProbe.cpp ${HIT_TEST_BOUNDARY_SOURCES})
        ADD_NEKTAR_EXECUTABLE(HITFieldCompare COMPONENT tests DEPENDS SolverUtils
            SOURCES ${HIT_TEST_DIR}/HITFieldCompare.cpp)
        ADD_NEKTAR_EXECUTABLE(HITFlowDiagnostics COMPONENT tests DEPENDS SolverUtils
            SOURCES ${HIT_TEST_DIR}/HITFlowDiagnostics.cpp)
        ADD_NEKTAR_EXECUTABLE(HITInitialField COMPONENT tests DEPENDS SolverUtils
            SOURCES ${HIT_TEST_DIR}/HITInitialField.cpp ${HIT_TEST_BOUNDARY_SOURCES})
        ADD_NEKTAR_EXECUTABLE(HITInflowBenchmark COMPONENT tests DEPENDS SolverUtils
            SOURCES ${HIT_TEST_DIR}/HITInflowBenchmark.cpp ${HIT_TEST_BOUNDARY_SOURCES})
        FOREACH(hit_tool HITBoundaryProbe HITFieldCompare HITFlowDiagnostics
                         HITInitialField HITInflowBenchmark)
            TARGET_INCLUDE_DIRECTORIES(${hit_tool} PRIVATE ${CMAKE_SOURCE_DIR}/solvers)
            TARGET_COMPILE_DEFINITIONS(${hit_tool} PRIVATE ENABLE_NEKTAR_EXCEPTIONS)
        ENDFOREACH()
        ADD_NEKTAR_TEST(HITInflow_m5_single_mode)
        ADD_NEKTAR_TEST(HITInflowSVVFourierPair)
        ADD_NEKTAR_TEST(HITInflowFieldCompareAbsoluteNorm)
        ADD_NEKTAR_TEST(HITInflowSeededH1D)
        ADD_NEKTAR_TEST(HITInflowTabulatedH1D)
        ADD_NEKTAR_TEST(HITInflowRestartH1D)
        ADD_NEKTAR_TEST(HITInflowProjectionH1D)
        ADD_NEKTAR_TEST(HITInflowDivergenceFlux3D)
        IF(NEKTAR_USE_FFTW)
            ADD_NEKTAR_TEST(HITInflowProjectionH1D_FFTW)
        ENDIF()
        IF(NEKTAR_USE_MPI)
            ADD_NEKTAR_TEST(HITInflowInitialFieldMPI)
            SET_TESTS_PROPERTIES(IncNavierStokesSolver_HITInflowInitialFieldMPI
                PROPERTIES PROCESSORS 4)
            ADD_NEKTAR_TEST(HITInflowMPIH1D)
            SET_TESTS_PROPERTIES(IncNavierStokesSolver_HITInflowMPIH1D
                PROPERTIES PROCESSORS 8)
        ENDIF()
    ENDIF()
    ADD_NEKTAR_TEST(ThermStab_dir)
    ADD_NEKTAR_TEST(ChannelStability)
    ADD_NEKTAR_TEST(ThermStab_adj)
    ADD_NEKTAR_TEST(ChanStabilityMaskFunc)
    ADD_NEKTAR_TEST(Channel_FlowRate_nonzeroBCs)
    ADD_NEKTAR_TEST(ChannelStabilityAdjoint) 
    # ADD_NEKTAR_TEST(ChanStability_Coupled_3D)
    ADD_NEKTAR_TEST(2DFlow_lineforcing_bcfromfile)
    ADD_NEKTAR_TEST(BercovierEngelman)
    ADD_NEKTAR_TEST(ChanFlow2D_bcsfromfiles)
    ADD_NEKTAR_TEST(ChanFlow3D_infTurb)
    ADD_NEKTAR_TEST(ChanFlow_3DH1D_bcsfromfiles)
    ADD_NEKTAR_TEST(ChannelFlowQuasi3D LENGTHY)
    ADD_NEKTAR_TEST(ChanFlow_3DH2D_MVM LENGTHY)
    ADD_NEKTAR_TEST(ChanFlow_3DH1D_pRef)
    ADD_NEKTAR_TEST(ChanFlow_LinNS_m8)
    ADD_NEKTAR_TEST(ChannelFlow2D)
    ADD_NEKTAR_TEST(ChanFlow_m3_ConOBC)
    ADD_NEKTAR_TEST(ChanFlow_m3_VCSWeakPress)
    ADD_NEKTAR_TEST(ChanFlow_m3_VCSWeakPress_ConOBC)
    ADD_NEKTAR_TEST(ChanFlow_m3_ConOBC_VCSImplicit)
    ADD_NEKTAR_TEST(ChanFlow_m3_SKS)
    ADD_NEKTAR_TEST(ChanFlow_m8)
    ADD_NEKTAR_TEST(ChanFlow_Accelerate)
    ADD_NEKTAR_TEST(ChanFlow_m8_BodyForce)
    ADD_NEKTAR_TEST(ChanFlow_m8_singular)
    ADD_NEKTAR_TEST(ChanFlow_V8P7_Avg)
    ADD_NEKTAR_TEST(ChanFlow_Tri_Lagrange)
    ADD_NEKTAR_TEST(Channel_Flow_3modes_rad)
    ADD_NEKTAR_TEST(channelTemp)
    ADD_NEKTAR_TEST(channelFlow_Heating_FlowRate)
    ADD_NEKTAR_TEST(channelFlow_nonuniform_Heating)
    ADD_NEKTAR_TEST(Couette_3DH2D_MVM)
    ADD_NEKTAR_TEST(FreeFallCyl)
    ADD_NEKTAR_TEST(Hex_channel_m3)
    ADD_NEKTAR_TEST(Hex_channel_varP)
    ADD_NEKTAR_TEST(Pyr_channel_m3)
    ADD_NEKTAR_TEST(Pyr_channel_varP)
    ADD_NEKTAR_TEST(Hex_channel_m6_nodalRestart)
    ADD_NEKTAR_TEST(Hex_channel_m3_srhs LENGTHY)
    ADD_NEKTAR_TEST(Hex_Kovasnay_SubStep)
    ADD_NEKTAR_TEST(PrismHex_channel_m4)
    ADD_NEKTAR_TEST(HM)
    ADD_NEKTAR_TEST(HM_Adj)
    ADD_NEKTAR_TEST(HM_Adj_SingleMode_files)
    ADD_NEKTAR_TEST(bfs_tg LENGTHY)
    ADD_NEKTAR_TEST(SM)
    ADD_NEKTAR_TEST(SM_base_file)
    ADD_NEKTAR_TEST(SM_Adj)
    ADD_NEKTAR_TEST(SM_Adj_base_file)
    ADD_NEKTAR_TEST(KovaFlow_m3)
    ADD_NEKTAR_TEST(KovaFlow_m3_GJP)
    ADD_NEKTAR_TEST(KovaFlow_m3_GJP_VCSWeakPressure)
    ADD_NEKTAR_TEST(KovaFlow_m3_GJP_VCSImplicit)
    ADD_NEKTAR_TEST(KovaFlow_m3_gmresZeroInputCheck)
    ADD_NEKTAR_TEST(KovaFlow_Tet_Curved_GJP_SemiImplicit)
    ADD_NEKTAR_TEST(KovaFlow_Tet_Curved_GJP_Explicit)
    ADD_NEKTAR_TEST(KovaFlow_m8)
    ADD_NEKTAR_TEST(KovaFlow_m8_short_ConOBC)
    ADD_NEKTAR_TEST(KovaFlow_m8_short_ConOBC_VCSWeakPress)
    ADD_NEKTAR_TEST(KovaFlow_m8_short_ConOBC_VCSImplicit)
    ADD_NEKTAR_TEST(KovaFlow_m8_short_HOBC)
    ADD_NEKTAR_TEST(KovaFlow_m8_short_HOBC_VCSWeakPress)
    ADD_NEKTAR_TEST(KovaFlow_m8_short_HOBC_VCSImplicit)
    ADD_NEKTAR_TEST(KovaFlow_m10_VCSImplicit)
    ADD_NEKTAR_TEST(KovaFlow_m10_VCSImplicit_MatFree)
    ADD_NEKTAR_TEST(KovaFlow_m10_taylorHood_VCSImplicit)
    ADD_NEKTAR_TEST(KovaFlow_m10_taylorHood_VCSImplicitFlexible)
    ADD_NEKTAR_TEST(KovaFlow_m10_taylorHood_VCSImplicitLoc)
    ADD_NEKTAR_TEST(KovaFlow_m10_taylorHood_VCSImplicitFlexibleLoc)
    ADD_NEKTAR_TEST(KovaFlow_m10_QuadTri_VCSImplicit)
    ADD_NEKTAR_TEST(KovaFlow_m10_QuadTri_VCSImplicit_Updated)
    ADD_NEKTAR_TEST(KovaFlow_m10_VCSImplicit_SVV)
    ADD_NEKTAR_TEST(KovaFlow_m10_3DH1D_absorption)
    ADD_NEKTAR_TEST(KovaFlow_varP)
    ADD_NEKTAR_TEST(KovaFlow_varP_per)
    #ADD_NEKTAR_TEST(KovaFlow_Oseen_m8)
    ADD_NEKTAR_TEST(KovaFlow_3DH1D_P5_20modes_MVM LENGTHY)
    ADD_NEKTAR_TEST(KovaFlow_3DH1D_P5_20modes_MVM_SVVHomo1D LENGTHY)
    ADD_NEKTAR_TEST(KovaFlow_3DH1D_P5_20modes_MVM_Deal LENGTHY)
    ADD_NEKTAR_TEST(KovaFlow_3DH1D_P5_20modes_SKS_MVM LENGTHY)
    ADD_NEKTAR_TEST(KovaFlow_3DH1D_P8_8modes_FFTW_VCSWeakPressure)
    ADD_NEKTAR_TEST(KovaFlow_SubStep_2order)
    ADD_NEKTAR_TEST(KovaFlow_SubStep_2order_VCSWeakPress)
    ADD_NEKTAR_TEST(Kovasznay3D_AllElmts_GJP)
    ADD_NEKTAR_TEST(Kovas_Quad6_Tri4_mixedbcs)
    ADD_NEKTAR_TEST(MovingRefFrame_translation)
    ADD_NEKTAR_TEST(MovingRefFrame_Rot_SimpleDomain)
    ADD_NEKTAR_TEST(MovingRefFrame_Rot_naca0012)
    ADD_NEKTAR_TEST(SelfPropel_MRF_naca0012)
    ADD_NEKTAR_TEST(Rayleigh_Bernard)
    ADD_NEKTAR_TEST(SinCos_LinNS_3DHom1D)
    ADD_NEKTAR_TEST(TaylorVor_dt1)
    ADD_NEKTAR_TEST(TaylorVor_dt2 LENGTHY)
    ADD_NEKTAR_TEST(Tet_Kovasnay_SubStep)
    ADD_NEKTAR_TEST(Tet_Kovasnay_HOBC)
    ADD_NEKTAR_TEST(Tet_Kovasnay_SVV_DGKer)
    ADD_NEKTAR_TEST(Tet_Kovasnay_PreconditionerBlock)
    ADD_NEKTAR_TEST(Tet_Kovasnay_PreconditionerLowEnergyBlock)
    ADD_NEKTAR_TEST(Tet_Kovasnay_GJP)
    ADD_NEKTAR_TEST(Tet_Kovasnay_GJP_VCSImplicit)
    ADD_NEKTAR_TEST(Tet_channel_m3)
    ADD_NEKTAR_TEST(Tet_channel_m3_ConOBC)
    ADD_NEKTAR_TEST(Tet_channel_m3_VCSWeakPress)
    ADD_NEKTAR_TEST(Tet_channel_m3_VCSWeakPress_ConOBC)
    ADD_NEKTAR_TEST(Tet_channel_m3_dealiasing)
    ADD_NEKTAR_TEST(Tet_channel_robin_m3)
    ADD_NEKTAR_TEST(Tet_channel_m4_per)
    ADD_NEKTAR_TEST(Tet_channel_varP)
    ADD_NEKTAR_TEST(TriQuadChannel)
    ADD_NEKTAR_TEST(Tet_prism_channel_rad LENGTHY)
    ADD_NEKTAR_TEST(Tet_equitri LENGTHY)
    ADD_NEKTAR_TEST(Prism_channel_m6 LENGTHY)
    ADD_NEKTAR_TEST(Prism_channel_m6_dealiasing LENGTHY)
    ADD_NEKTAR_TEST(Prism_channel_varP)
    ADD_NEKTAR_TEST(ChannelSpongeNSE)
    ADD_NEKTAR_TEST(ChannelSpongeLNSE)
    ADD_NEKTAR_TEST(ChanFlow_Standard_BodyForce)
    ADD_NEKTAR_TEST(Cyl_AdaptiveSFD)
    ADD_NEKTAR_TEST(Womersley_PipeFlow)
    ADD_NEKTAR_TEST(CylFlow_adaptiveP LENGTHY)
    ADD_NEKTAR_TEST(CylFlow_Mov_mapping)
    ADD_NEKTAR_TEST(PPF_R10000_ModifiedArnoldi_Shift)
    ADD_NEKTAR_TEST(PPF_R15000_ModifiedArnoldi_Shift)
    ADD_NEKTAR_TEST(ChanFlow_m8_Flowrate)
    ADD_NEKTAR_TEST(ChanFlow_3DH1D_FlowrateExplicit_MVM)
    ADD_NEKTAR_TEST(PhaseSampling)
    ADD_NEKTAR_TEST(CylFlow2D_SPM)
    ADD_NEKTAR_TEST(CylFlow2DMov_SPM)
    ADD_NEKTAR_TEST(CylFlow3DH1D_SPM)
    ADD_NEKTAR_TEST(manufactured_3DH1D)
    ADD_NEKTAR_TEST(PlungingAirfoil_2D)
    ADD_NEKTAR_TEST(3DH1D_L2ErrorFilter)

    IF(NOT WIN32 AND NEKTAR_UTILITY_FIELDCONVERT)
        ADD_NEKTAR_TEST(ScaleReynoldsShearFilter_ThermalChannel2D)
        ADD_NEKTAR_TEST(ScaleReynoldsShearFilter_ThermalChannel3DH1D)
    ENDIF()

    IF(NEKTAR_USE_SCOTCH)
        # This tests exhibits sensitivity to the choice of DirectStaticCond vs.
        # DirectMultiLevelStaticCond.
        ADD_NEKTAR_TEST(Pyr_channel_SVV)
    ENDIF()

    IF(NEKTAR_USE_ARPACK)
        ADD_NEKTAR_TEST(ChanStability_adj_Ar)
        ADD_NEKTAR_TEST(ThermStab_dir_Ar)
        ADD_NEKTAR_TEST(ThermStab_adj_Ar)
        ADD_NEKTAR_TEST(bfs_tg-AR)
        ADD_NEKTAR_TEST(ChanStability_Coupled)
        ADD_NEKTAR_TEST(PPF_R10000_Arpack_LM)
        ADD_NEKTAR_TEST(PPF_R15000_Arpack_Shift)
        # ADD_NEKTAR_TEST(PPF_R15000_Arpack_NoImagShift)
        ADD_NEKTAR_TEST(PPF_R15000_Arpack_NoImagShift_LM)
    ENDIF(NEKTAR_USE_ARPACK)

    IF(NEKTAR_USE_FFTW)
        ADD_NEKTAR_TEST(ChanFlow_3DH1D_FFT)
        ADD_NEKTAR_TEST(ChanFlow_3DH1D_FFT_ConOBC)
        ADD_NEKTAR_TEST(ChanFlow_3DH2D_FFT)
        ADD_NEKTAR_TEST(CylFlow_MovBody)
        ADD_NEKTAR_TEST(KovaFlow_3DH1D_MVM_FFTW_Consistency)
        ADD_NEKTAR_TEST(KovaFlow_3DH1D_P5_6modes_FFTW_MixedDeal)
        ADD_NEKTAR_TEST(KovaFlow_3DH1D_P8_16modes_Mapping-explicit)
        ADD_NEKTAR_TEST(KovaFlow_3DH1D_P8_16modes_Mapping-implicit)
        ADD_NEKTAR_TEST(KovaFlow_3DH2D_xy_FFT)
        ADD_NEKTAR_TEST(KovaFlow_3DH2D_xz_FFT)
        ADD_NEKTAR_TEST(KovaFlow_3DH2D_xy_FFT_MovRefFrame)
        ADD_NEKTAR_TEST(KovaFlow_m8_short_HOBC_3D1H)
        ADD_NEKTAR_TEST(KovaFlow_3DH1D_adaptive_16modes_FFTW_Mapping LENGTHY)
        ADD_NEKTAR_TEST(Couette_3DH2D_FFT)
    ENDIF(NEKTAR_USE_FFTW)

    IF(NEKTAR_USE_MPI)
        ADD_NEKTAR_TEST(ChanFlow_3DH1D_Parallel_mode1)
        ADD_NEKTAR_TEST(ChanFlow_3DH1D_Parallel_mode2)
        ADD_NEKTAR_TEST(ChanFlow_m3_par)
        ADD_NEKTAR_TEST(FlapAirfoil_3DH1D_par LENGTHY)
        ADD_NEKTAR_TEST(PlungingAirfoil_3DH1D_parBCs LENGTHY)
        ADD_NEKTAR_TEST(ChanFlow_m8_BodyForce_par LENGTHY)
        ADD_NEKTAR_TEST(Cyl_AdaptiveSFD_par)

        IF(NOT WIN32)
            ADD_NEKTAR_TEST(KovaFlow_expFromFile_par)
            ADD_NEKTAR_TEST(KovaFlow_varP_per_par)
            ADD_NEKTAR_TEST(Hex_channel_varP_par)
        ENDIF()

        ADD_NEKTAR_TEST(Hex_channel_m3_par LENGTHY)
        ADD_NEKTAR_TEST(Pyr_channel_m6_par LENGTHY)
        ADD_NEKTAR_TEST(Hex_channel_m3_srhs_par LENGTHY)

        IF(NOT WIN32)
            ADD_NEKTAR_TEST(Tet_channel_m4_per_xxt_ml_par)
        ENDIF()

        ADD_NEKTAR_TEST(Tet_channel_m8_par)
        ADD_NEKTAR_TEST(Tet_channel_m8_iter_ml_par LENGTHY)
        ADD_NEKTAR_TEST(bfs_tg_par LENGTHY)

        ADD_NEKTAR_TEST(ChanFlow_m8_Flowrate_par)
        ADD_NEKTAR_TEST(ChanFlow_3DH1D_Flowrate_MVM)
        ADD_NEKTAR_TEST(ChanFlow_3DH1D_Flowrate_MVM_hybrid)
        ADD_NEKTAR_TEST(Hex_channel_m3_Flowrate)
        ADD_NEKTAR_TEST(ChanFlow_3DH1D_FlowrateExplicit_MVM_par)
        ADD_NEKTAR_TEST(ChanFlow_3DH1D_FlowrateExplicit_MVM_par_hybrid)
        ADD_NEKTAR_TEST(PointsTracking2D)
        ADD_NEKTAR_TEST(PointsTracking3D)
        ADD_NEKTAR_TEST(PointsTracking3DH1D)

        IF(NEKTAR_USE_PETSC)
            ADD_NEKTAR_TEST(Tet_channel_m8_petsc_sc_par LENGTHY)
        ENDIF(NEKTAR_USE_PETSC)
        IF(NEKTAR_USE_HDF5)
            ADD_NEKTAR_TEST(TriQuadChannelHdf5)
            ADD_NEKTAR_TEST(Tet_channel_m4_per_hdf5)
            ADD_NEKTAR_TEST(Tet_channel_m8_hdf)
            ADD_NEKTAR_TEST(KovaFlow_3DH1D_P5_20modes_MVM_hdf5)
            ADD_NEKTAR_TEST(CylFlow2D_FieldConvertFilter_Hdf5)
            ADD_NEKTAR_TEST(CylFlow2D_CheckpointFilter_Hdf5)
            ADD_NEKTAR_TEST(Tet_hdf5_multilevel_partition)
            IF(NEKTAR_USE_SCOTCH)
                ADD_NEKTAR_TEST(ChannelFlow3D)
            ENDIF()
        ENDIF(NEKTAR_USE_HDF5)
    ENDIF(NEKTAR_USE_MPI)
ENDIF(NEKTAR_SOLVER_INCNAVIERSTOKES)
