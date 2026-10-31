///////////////////////////////////////////////////////////////////////////////
// File: HITInflow.h
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_INCNAVIERSTOKES_HITINFLOW_H
#define NEKTAR_INCNAVIERSTOKES_HITINFLOW_H

#include <IncNavierStokesSolver/BoundaryConditions/HITInflowConfig.h>
#include <IncNavierStokesSolver/BoundaryConditions/IncBaseCondition.h>
#include <SolverUtils/Turbulence/HITPlaneEvaluator.h>

namespace Nektar
{
/// Three-component trace of a frozen, solenoidal Fourier realization.
class HITInflow : public IncBaseCondition
{
public:
    friend class MemoryManager<HITInflow>;

    static IncBaseConditionSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr session,
        Array<OneD, MultiRegions::ExpListSharedPtr> fields,
        Array<OneD, SpatialDomains::BoundaryConditionShPtr> conditions,
        Array<OneD, MultiRegions::ExpListSharedPtr> expansions, int boundary,
        int spaceDimension, int boundaryDimension)
    {
        auto result = MemoryManager<HITInflow>::AllocateSharedPtr(
            session, fields, conditions, expansions, boundary, spaceDimension,
            boundaryDimension);
        result->Initialise(session);
        return result;
    }

    static std::string className;

    /// The manager supplies the descriptor after all-rank validation. Neither
    /// this method nor the local constructor starts a session collective.
    void Configure(
        std::shared_ptr<const HITInflowConfiguration> configuration,
        std::shared_ptr<const SolverUtils::HITGenerator> generator);

protected:
    HITInflow(const LibUtilities::SessionReaderSharedPtr session,
              Array<OneD, MultiRegions::ExpListSharedPtr> fields,
              Array<OneD, SpatialDomains::BoundaryConditionShPtr> conditions,
              Array<OneD, MultiRegions::ExpListSharedPtr> expansions,
              int boundary, int spaceDimension, int boundaryDimension);

    void v_Initialise(
        const LibUtilities::SessionReaderSharedPtr &session) override;
    void v_UpdatePrescribedVelocity(NekDouble time) override;

private:
    std::shared_ptr<const HITInflowConfiguration> m_configuration;
    std::shared_ptr<const SolverUtils::HITGenerator> m_generator;
    std::unique_ptr<SolverUtils::HITPlaneEvaluator> m_evaluator;
    std::vector<SolverUtils::HITVector> m_velocities;
    int m_reportEvaluator = 0;
    int m_requireAccelerated = 0;
    int m_rank = 0;
};
} // namespace Nektar

#endif
