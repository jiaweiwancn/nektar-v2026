///////////////////////////////////////////////////////////////////////////////
// File: HITInflowGeometry.h
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_INCNAVIERSTOKES_HITINFLOWGEOMETRY_H
#define NEKTAR_INCNAVIERSTOKES_HITINFLOWGEOMETRY_H

#include <IncNavierStokesSolver/BoundaryConditions/HITInflowConfig.h>
#include <MultiRegions/ExpList.h>
#include <SpatialDomains/Conditions.h>
#include <vector>

namespace Nektar
{
class HITInflowGeometry
{
public:
    /// Collective on every session rank, including ranks without inlet faces.
    /// The supplied BoundaryConditions must already have been constructed on
    /// every rank. Return global XML REF IDs in local boundary-array order.
    /// Validate six explicit box sides in full 3D, or four xy sides and the
    /// implicit Fourier z period in Homogeneous1D.
    static std::vector<int> Validate(
        const LibUtilities::SessionReaderSharedPtr &session,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const HITInflowConfiguration &configuration,
        const SpatialDomains::BoundaryConditions &boundaryConditions);
};
} // namespace Nektar

#endif
