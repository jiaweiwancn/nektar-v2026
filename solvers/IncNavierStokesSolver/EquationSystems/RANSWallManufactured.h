// SPDX-License-Identifier: MIT
#ifndef NEKTAR_RANS_WALL_MANUFACTURED_H
#define NEKTAR_RANS_WALL_MANUFACTURED_H
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
namespace Nektar
{
bool RANSUsesAcceptedManufactured(const LibUtilities::SessionReaderSharedPtr &);
std::string RANSAcceptedManufacturedIdentity(const LibUtilities::SessionReaderSharedPtr &);
struct RANSWallManufacturedSample
{
    std::array<double,5> value{}; // Ux,Uz,k,total epsilon,augmented pressure
    std::array<SolverUtils::RANSVector,5> gradient{};
    SolverUtils::RANSVector body{},velocityTime{},normal{},tangent{};
    SolverUtils::RANSTensor reynolds{},quadratic{},stress{},velocityGradient{};
    double kSource=0,epsilonSource=0,kTime=0,epsilonTime=0;
    double kDiffusion=0,epsilonDiffusion=0,production=0,epsilonProduction=0,epsilonDestruction=0;
    double K=0,cutoff=0,J=0,a=0,gamma=0,radicand=0;
};
/** Continuous manufactured data only. Independent Taylor derivatives and
 * planar constitutive algebra; never consumes a native FE residual/state.
 * The physical production configuration cannot select this provider. */
class RANSWallManufactured
{
    struct Impl;
    std::shared_ptr<Impl> impl;
public:
    explicit RANSWallManufactured(const RANSWallEvolutionPatch &);
    RANSWallManufacturedSample At(const RANSWallCoordinate &,double,bool) const;
    RANSWallEvolutionData Data() const;
    std::vector<RANSWallEvolutionFunction> InitialFunctions() const;
};
}
#endif
