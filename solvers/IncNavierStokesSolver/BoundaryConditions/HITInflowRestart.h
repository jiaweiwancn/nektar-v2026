///////////////////////////////////////////////////////////////////////////////
// File: HITInflowRestart.h
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_INCNAVIERSTOKES_HITINFLOWRESTART_H
#define NEKTAR_INCNAVIERSTOKES_HITINFLOWRESTART_H

#include <IncNavierStokesSolver/BoundaryConditions/HITInflowConfig.h>
#include <LibUtilities/BasicUtils/FieldIO.h>

namespace Nektar
{
class HITInflowRestart
{
public:
    /// Required realization identity and numerical reconstruction metadata.
    static LibUtilities::FieldMetaDataMap Describe(
        const HITInflowConfiguration &configuration,
        const SolverUtils::HITGenerator &generator);
    /// Pure comparison, also used by the boundary unit tests.
    static void Verify(const LibUtilities::FieldMetaDataMap &stored,
                       const LibUtilities::FieldMetaDataMap &expected);
    /// Collective, after native restart import and before the first trace/dump.
    static void Prepare(
        const LibUtilities::SessionReaderSharedPtr &session,
        const HITInflowConfiguration &configuration,
        const SolverUtils::HITGenerator &generator,
        LibUtilities::FieldMetaDataMap &metadata, NekDouble initialTime);
};
} // namespace Nektar

#endif
