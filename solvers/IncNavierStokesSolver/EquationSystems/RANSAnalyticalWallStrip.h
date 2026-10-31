///////////////////////////////////////////////////////////////////////////////
// File: RANSAnalyticalWallStrip.h
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_RANS_ANALYTICAL_WALL_STRIP_H
#define NEKTAR_RANS_ANALYTICAL_WALL_STRIP_H

#include <MultiRegions/ExpList.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallBoundary.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSession.h>
#include <SolverUtils/Turbulence/AnalyticalWallFunction.h>
#include <SolverUtils/Turbulence/RANSWallStrip.h>
#include <fstream>

namespace Nektar
{
/** Flat CV-to-CG adapter. P is eliminated using the M1 matching operators;
 * no separate strip history is needed because these operators are linear and
 * are applied to the complete IMEX predictor as well as the current iterate.
 * Collective N sampling includes ranks with no local wall faces.
 */
class RANSAnalyticalWallStrip : public RANSWallBoundary
{
public:
    using Fields = Array<OneD, Array<OneD, NekDouble>>;
    using Gradients = Array<OneD, Fields>;
    RANSAnalyticalWallStrip(
        const LibUtilities::SessionReaderSharedPtr &session,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        int dimension, NekDouble wallOrigin,
        const SolverUtils::RANSNeutralReference &reference);
    void CheckInitialState(const Fields &state) const;
    std::string RestartSignature() const;
    void Update(const Fields &state, const Fields &predictor,
                const Gradients &gradient, const Fields &pressureForce,
                const Fields &convection, const Fields &bodyForce, NekDouble alpha);
    void AddRobinAuditResidual(int field, const Fields &state,
                              Array<OneD, NekDouble> &residual) const;
    void WritePointDiagnostics(NekDouble time, int phase);
    void WriteDiagnostics(NekDouble time, const Fields &state,
                          const Fields &gain, const Fields &loss,
                          const Array<OneD, NekDouble> &pressure);

private:
    Fields Trace(const Fields &values) const;
    Fields SampleN(const Fields &values) const;
    void SetBoundary(int field, const Array<OneD, NekDouble> &values,
                     bool dirichlet);
    void Require(bool valid, const std::string &message) const;
    NekDouble SurfaceIntegral(const Array<OneD, NekDouble> &values) const;
    NekDouble NaturalFlux(int field, bool &available) const;

    LibUtilities::SessionReaderSharedPtr m_session;
    Array<OneD, MultiRegions::ExpListSharedPtr> m_fields;
    SolverUtils::RANSWallStripGeometry m_geometry;
    SolverUtils::AWFGeometry m_awfGeometry;
    SolverUtils::AWFParameters m_parameters;
    SolverUtils::RANSNeutralReference m_reference;
    SolverUtils::AWFSwitchParameters m_switch;
    SolverUtils::AWFSpatialTransition m_spatial;
    Array<OneD,NekDouble> m_wallX, m_roughness, m_mask;
    SolverUtils::AWFFailurePolicy m_failurePolicy;
    std::string m_wallModel;
    bool m_equilibrium = false, m_switched = false, m_writePoints = false;
    NekDouble m_fixedWeight = 1, m_minimumRe = 100;
    struct PointData
    {
        SolverUtils::AWFIndependentSamples samples{};
        SolverUtils::EquilibriumAWFValues selected{};
        SolverUtils::EquilibriumAWFReference reference{};
        SolverUtils::AWFValues original{};
        double wallRe = 0;
    };
    std::vector<PointData> m_points;
    int m_dim, m_localPoints, m_globalPoints, m_globalOffset, m_remoteSamples;
    int m_wallRanks = 0, m_remoteWithoutWallSamples = 0;
    NekDouble m_minimumK, m_minimumEpsilon;
    NekDouble m_alpha = 0;
    std::vector<int> m_boundaries, m_offsets, m_sampleElements;
    std::vector<Array<OneD, NekDouble>> m_sampleCoordinates;
    Fields m_predictorP, m_g, m_pressureP, m_volumePredictor;
    Array<OneD, NekDouble> m_production, m_dissipation, m_epsilonP,
        m_epsilonH, m_kP, m_tractionX, m_tractionY, m_forcingX;
    std::ofstream m_diagnostics;
    std::ofstream m_pointDiagnostics;
};
} // namespace Nektar
#endif
