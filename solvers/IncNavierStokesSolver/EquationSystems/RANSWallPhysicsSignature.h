// SPDX-License-Identifier: MIT
#ifndef NEKTAR_RANS_WALL_PHYSICS_SIGNATURE_H
#define NEKTAR_RANS_WALL_PHYSICS_SIGNATURE_H
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <iomanip>
#include <sstream>
namespace Nektar
{
// Shared byte-for-byte native physical identity, also used for immutable
// checkpoint replays. Only the named projection resource parameters are omitted.
// The finite-QP numerical backend is a SolverInfo execution control, recorded
// in PimpleKBoundControls: it solves the same raw-centered mass QP. The floor,
// inventory, metric and mathematical correction policy retain strict identity.
// RANSWallDenseEvolutionCapacity/ResourceDesign are execution-only SolverInfo
// controls, reported separately. They do not enter patch.Signature or this
// Parameters loop; no additional physical parameter is excluded.
inline std::string RANSWallPhysicsSignature(const RANSWallEvolutionPatch &patch,
    const LibUtilities::SessionReaderSharedPtr &session,double eta,bool neutral,
    double ustar,const std::array<double,2> &direction,bool chartFunctions)
{
    std::ostringstream s;s<<std::setprecision(17)<<patch.Signature()<<":eta="<<eta<<":neutral="<<neutral<<":us="<<ustar<<":direction="<<direction[0]<<','<<direction[1]<<":chart_functions="<<chartFunctions;
    if(!neutral && patch.CappedDestruction() && !patch.Mapped() && !patch.TransitionEnabled())
        s<<":capped_initial_trace=1";
    if(patch.PimpleSingleSweep())s<<":temporal_coupling=PimpleSingleSweep"<<patch.PimpleKBoundSignature();
    for(const auto &[name,value]:session->GetParameters())
        if(!(patch.ColumnAverageK() && name=="RANSWALLSAMPLEHEIGHT") && !RANSWallEvolutionPatch::PimpleKBoundWorkParameter(name) && name!="NUMSTEPS" && name!="IO_CHECKSTEPS" && name!="IO_INFOSTEPS" && name!="IO_TIMER_LEVEL" && name!="RANSWALLOWNERSHIPSHIFT" && name!="RANSWALLPIMPLEREJECTSTEP" && name!="RANSWALLPIMPLEHARDAUDITREJECTSTEP" && name!="RANSWALLPIMPLEKFLOORRELATIVE" && name!="RANSWALLPIMPLEBUDGETTIME")s<<'|'<<name<<'='<<value;
    const std::vector<std::string> names{"InitialConditionsInner","InitialConditionsOuter","InitialPerturbationInner","InitialPerturbationOuter","BoundaryValuesInner","BoundaryValuesOuter",
        "BodyForceInner","BodyForceOuter","ScalarSourceInner","ScalarSourceOuter","ExactSolutionInner","ExactSolutionOuter",
        "OutletTractionInner","OutletTractionOuter","TopTractionOuter","TopDiffusiveFluxOuter"};
    for(const auto &name:names)for(const auto &variable:session->GetVariables())if(session->DefinesFunction(name,variable))
        s<<'|'<<name<<'/'<<variable<<'='<<session->GetFunction(name,variable)->GetExpression();
    return s.str();

}
}
#endif
