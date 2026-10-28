///////////////////////////////////////////////////////////////////////////////
// File: HITInflow.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/BoundaryConditions/HITInflow.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace Nektar
{
std::string HITInflow::className = GetIncBCFactory().RegisterCreatorFunction(
    "HITInflow", HITInflow::create,
    "Frozen homogeneous isotropic turbulence velocity inlet");

HITInflow::HITInflow(
    const LibUtilities::SessionReaderSharedPtr session,
    Array<OneD, MultiRegions::ExpListSharedPtr> fields,
    Array<OneD, SpatialDomains::BoundaryConditionShPtr> conditions,
    Array<OneD, MultiRegions::ExpListSharedPtr> expansions, int boundary,
    int spaceDimension, int boundaryDimension)
    : IncBaseCondition(session, fields, conditions, expansions, boundary,
                       spaceDimension, boundaryDimension)
{
    for (int component = 0; component < 3; ++component)
    {
        m_BndConds[component] = conditions[component];
        m_BndExp[component]   = expansions[component];
    }
}

void HITInflow::v_Initialise(
    const LibUtilities::SessionReaderSharedPtr &session)
{
    IncBaseCondition::v_Initialise(session);
    m_coords = Array<OneD, Array<OneD, NekDouble>>(3);
    for (int component = 0; component < 3; ++component)
    {
        m_coords[component] = Array<OneD, NekDouble>(m_npoints);
    }
    m_BndExp[0]->GetCoords(m_coords[0], m_coords[1], m_coords[2]);
}

void HITInflow::Configure(
    std::shared_ptr<const HITInflowConfiguration> configuration,
    std::shared_ptr<const SolverUtils::HITGenerator> generator)
{
    if (!configuration || !generator || m_spacedim != 3 || m_bnddim != 3)
    {
        throw std::invalid_argument(
            "HITInflow requires a configured full three-dimensional inlet");
    }
    for (int component = 0; component < 3; ++component)
    {
        if (m_BndExp[component]->GetNpoints() != m_npoints ||
            m_BndExp[component]->GetNcoeffs() != m_BndExp[0]->GetNcoeffs())
        {
            throw std::invalid_argument(
                "HITInflow velocity components must use matching expansions");
        }
        Array<OneD, Array<OneD, NekDouble>> coordinates(3);
        for (int axis = 0; axis < 3; ++axis)
        {
            coordinates[axis] = Array<OneD, NekDouble>(m_npoints);
        }
        m_BndExp[component]->GetCoords(coordinates[0], coordinates[1],
                                       coordinates[2]);
        for (int axis = 0; axis < 3; ++axis)
        {
            for (int point = 0; point < m_npoints; ++point)
            {
                if (coordinates[axis][point] != m_coords[axis][point])
                {
                    throw std::invalid_argument(
                        "HITInflow velocity quadrature coordinates must match");
                }
            }
        }
    }
    m_configuration = std::move(configuration);
    m_generator     = std::move(generator);
}

void HITInflow::v_UpdatePrescribedVelocity(NekDouble time)
{
    ASSERTL0(m_generator, "HITInflow was not configured before its first update");
    for (int point = 0; point < m_npoints; ++point)
    {
        const SolverUtils::HITVector location{
            m_coords[0][point], m_coords[1][point], m_coords[2][point]};
        const auto sample = m_generator->Evaluate(location, time);
        for (int component = 0; component < 3; ++component)
        {
            m_BndExp[component]->UpdatePhys()[point] = sample.velocity[component];
        }
    }
    for (int component = 0; component < 3; ++component)
    {
        auto &boundary = m_BndExp[component];
        boundary->SetWaveSpace(false);
        boundary->FwdTransBndConstrained(boundary->GetPhys(),
                                         boundary->UpdateCoeffs());
    }
}
} // namespace Nektar
