# Shared by the normal production target and isolated standard-main build.
# The identity covers numerical sources and headers, independent of output paths.
function(wj_h27_source_identity root output)
    file(GLOB wj_native_files "${root}/solvers/IncNavierStokesSolver/EquationSystems/WJNative*.h"
                              "${root}/solvers/IncNavierStokesSolver/EquationSystems/WJNative*.cpp")
    list(APPEND wj_native_files
        "${root}/library/SolverUtils/Turbulence/WJEARSM.h"
        "${root}/library/SolverUtils/Turbulence/WJEARSM.cpp"
        "${root}/library/SolverUtils/Turbulence/RANSModels.h"
        "${root}/library/SolverUtils/Turbulence/RANSModels.cpp"
        "${root}/library/SolverUtils/Turbulence/AnalyticalWallFunction.h"
        "${root}/library/SolverUtils/Turbulence/AnalyticalWallFunction.cpp")
    list(SORT wj_native_files)
    set(binding "WJH27-v2|C++17|no-fast-math|fp-contract=off")
    foreach(path IN LISTS wj_native_files)
        file(SHA256 "${path}" digest)
        string(APPEND binding "|${digest}")
    endforeach()
    string(SHA256 identity "${binding}")
    set(${output} "${identity}" PARENT_SCOPE)
endfunction()
