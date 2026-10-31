///////////////////////////////////////////////////////////////////////////////
// File: HITInflow.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/BoundaryConditions/HITInflow.h>
#include <algorithm>
#include <cmath>
#include <iostream>
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
    session->LoadParameter("HITReportEvaluator", m_reportEvaluator, 0);
    session->LoadParameter("HITRequireAccelerated", m_requireAccelerated, 0);
    m_rank = session->GetComm()->GetRank();
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
    const bool supported = m_field->GetExpType() == MultiRegions::e3D ||
                           m_field->GetExpType() == MultiRegions::e3DH1D;
    if (!configuration || !generator || !supported || m_spacedim != 3 ||
        (m_bnddim != 2 && m_bnddim != 3))
    {
        throw std::invalid_argument(
            "HITInflow requires a configured full 3D or Homogeneous1D inlet");
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
    std::vector<SolverUtils::HITVector> points(m_npoints);
    for (int point = 0; point < m_npoints; ++point)
    {
        points[point] = {m_coords[0][point], m_coords[1][point],
                         m_coords[2][point]};
    }
    m_evaluator = std::make_unique<SolverUtils::HITPlaneEvaluator>(
        *m_generator, points);
    ASSERTL0(!m_requireAccelerated || m_evaluator->IsGrouped(),
             "HITRequireAccelerated: bounded plane acceleration is unavailable");
    if (m_reportEvaluator)
    {
        std::cout << "HIT_EVALUATOR rank " << m_rank << " points " << m_npoints
                  << " method " << m_evaluator->Method() << " cache_bytes "
                  << m_evaluator->CacheBytes() << " modes "
                  << m_generator->GetModes().size() << std::endl;
    }
}

void HITInflow::v_UpdatePrescribedVelocity(NekDouble time)
{
    ASSERTL0(m_evaluator, "HITInflow was not configured before its first update");
    m_evaluator->Evaluate(time, m_velocities);
    for (int point = 0; point < m_npoints; ++point)
    {
        for (int component = 0; component < 3; ++component)
        {
            m_BndExp[component]->UpdatePhys()[point] = m_velocities[point][component];
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
