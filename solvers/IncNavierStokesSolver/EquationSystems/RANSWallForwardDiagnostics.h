// SPDX-License-Identifier: MIT
#ifndef NEKTAR_RANS_WALL_FORWARD_DIAGNOSTICS_H
#define NEKTAR_RANS_WALL_FORWARD_DIAGNOSTICS_H
#include <SolverUtils/Turbulence/RANSModels.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>

namespace Nektar::RANSWallDiagnostics
{
struct ForwardObservation
{
    SolverUtils::RANSClosureValues closure;
    std::array<double,3> covarianceEigenvalues;
    bool tangentAvailable=false;
    const char *tangentStatus="unavailable_forward_only_new_model";
};
// This is the production forward-audit path, also exposed to the deterministic
// audit regression. A diagnostic stencil is not part of physical admissibility.
inline ForwardObservation Observe(const SolverUtils::RANSResolvedModel &model,
    double k,double epsilon,const SolverUtils::RANSTensor &gradient,double distance)
{
    if(!model.NewModel())throw std::invalid_argument("new-model forward audit required");
    ForwardObservation result;
    result.closure=SolverUtils::EvaluateRANSClosure(model,k,epsilon,gradient,distance);
    result.covarianceEigenvalues=RANSWallSolveJson::Eigenvalues(result.closure.reynoldsStress);
    return result;
}
}
#endif
