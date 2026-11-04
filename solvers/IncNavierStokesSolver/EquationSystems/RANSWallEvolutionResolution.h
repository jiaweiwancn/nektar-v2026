///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_RANS_WALL_EVOLUTION_RESOLUTION_H
#define NEKTAR_RANS_WALL_EVOLUTION_RESOLUTION_H
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>

namespace Nektar
{
/** Physical spatial terms. Auxiliary E/C and alternative arithmetic entries
 * are deliberately absent: destruction is assembled once as +D in -Sk. */
const std::vector<std::string> &RANSWallResolutionSpatialTerms();

struct RANSWallScalarMassRate
{
    RANSWallAlgebra::Vector rate,boundaryRate,spatialLoad,freeRows,fixedRows,
        freeMassResidual;
    std::map<std::string,RANSWallAlgebra::Vector> componentRates;
    double contactValue=0,contactRate=0;
    std::map<std::string,double> contactComponents;
};
/** Consistent native physical mass, never lumped. The boundary rate is an
 * explicit input; a missing vector means the declared static boundary rate0.
 * This k-only instantaneous equation does not itself solve the dependent EPS
 * DAE. That separate coupling/solvability obligation is retained by callers. */
RANSWallScalarMassRate RANSWallResolutionMassRate(const RANSWallEvolutionPatch &patch,
    const RANSWallEvolutionPatch::ResidualResolutionResult &resolution,int block,
    const RANSWallCoordinate &contact,const RANSWallAlgebra::Vector &boundaryRate={});
}
#endif
