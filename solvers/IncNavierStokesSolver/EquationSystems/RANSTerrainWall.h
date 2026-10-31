///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_RANS_TERRAIN_WALL_H
#define NEKTAR_RANS_TERRAIN_WALL_H
#include <IncNavierStokesSolver/EquationSystems/RANSWallBoundary.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSession.h>
#include <SolverUtils/Turbulence/RANSWallStrip.h>
#include <fstream>

namespace Nektar
{
/** Normal-coordinate, thin-layer CV/CG coupling for a 2D terrain graph (and
 * its spanwise extrusion). Cartesian vector Robin rows solve cross-component
 * terms together. A boundary multiplier imposes the normal
 * constraint in the continuous free velocity trace space; augmentation affects
 * convergence only. All boundary loads are per actual CG interface area.
 */
class RANSTerrainWall : public RANSWallBoundary
{
public:
    RANSTerrainWall(const LibUtilities::SessionReaderSharedPtr &session,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields, int dimension,
        const SolverUtils::RANSWallGeometry &geometry,
        const SolverUtils::RANSNeutralReference &reference);
    void CheckInitialState(const Fields &) const override;
    std::string RestartSignature() const override;
    void Update(const Fields &, const Fields &, const Gradients &,
                const Fields &, const Fields &, const Fields &, NekDouble) override;
    void AddRobinAuditResidual(int, const Fields &, Array<OneD, NekDouble> &) const override;
    void WritePointDiagnostics(NekDouble, int) override;
    void WriteDiagnostics(NekDouble, const Fields &, const Fields &,
                          const Fields &, const Array<OneD, NekDouble> &) override;
    NekDouble ConstraintResidual() const override { return m_constraint; }
    Array<OneD,NekDouble> IterationState() const override { return m_multiplier; }
    void SetIterationState(const Array<OneD,NekDouble> &) override;
    Array<OneD,NekDouble> IterationMetric() const override;
    NekDouble IterationReactionScale() const override { return m_beta*std::sqrt(m_normalRelaxation); }
    void SetTrialMode(bool value) override { m_trialMode = value; }

private:
    void UpdateContinuity(const Fields &,const Fields &);
    Fields Trace(const Fields &) const;
    Fields SampleN(const Fields &) const;
    void SetBoundary(int, const Array<OneD, NekDouble> &, bool);
    void Require(bool, const std::string &) const;
    Array<OneD, NekDouble> ProjectNormal(const Array<OneD, NekDouble> &) const;
    Array<OneD, NekDouble> NormalLoad(const Array<OneD, NekDouble> &) const;
    Array<OneD, NekDouble> NormalValues(const Array<OneD, NekDouble> &) const;
    void BuildNormalDirect();
    Array<OneD,NekDouble> ProjectConstraint(const Array<OneD,NekDouble> &) const;
    Array<OneD,NekDouble> AdjointConstraint(const Array<OneD,NekDouble> &) const;
    NekDouble SurfaceIntegral(const Array<OneD, NekDouble> &) const;
    SolverUtils::RANSVector Vector(const Fields &, int) const;
    SolverUtils::AWFVector Tangent(const SolverUtils::RANSVector &, int) const;
    SolverUtils::RANSVector Cartesian(const SolverUtils::AWFVector &, int) const;
    NekDouble Component(const SolverUtils::RANSVector &, int) const;

    LibUtilities::SessionReaderSharedPtr m_session;
    Array<OneD, MultiRegions::ExpListSharedPtr> m_fields;
    int m_dim, m_localPoints = 0, m_globalPoints = 0, m_globalOffset = 0;
    int m_remoteSamples = 0;
    NekDouble m_samplingResidual = 0;
    NekDouble m_projectionAudit = 0;
    NekDouble m_resolvedVolume = 0, m_unresolvedVolume = 0;
    std::string m_wallModel, m_tag;
    bool m_equilibrium = false, m_log = false, m_switched = false;
    SolverUtils::RANSWallGeometry m_terrain;
    SolverUtils::RANSWallStripGeometry m_strip;
    SolverUtils::AWFGeometry m_awf;
    SolverUtils::AWFParameters m_parameters;
    SolverUtils::RANSNeutralReference m_reference;
    SolverUtils::AWFSpatialTransition m_spatial;
    SolverUtils::AWFSwitchParameters m_switch;
    SolverUtils::AWFFailurePolicy m_failurePolicy;
    double m_fixedWeight = 1, m_minimumRe = 100, m_augmentation = 2, m_normalRelaxation = .5;
    double m_resolvedNu = 0;
    bool m_multiplierInitialized = false;
    bool m_trialMode = false;
    double m_maximumCurvature = .25, m_geometryTolerance = 1e-7;
    double m_z0 = .02, m_shift = 0, m_minimumK = 0, m_minimumEpsilon = 0;
    NekDouble m_constraint = 0;
    std::vector<int> m_boundaries, m_offsets, m_sampleElements;
    std::vector<int> m_traceOffsets;
    Array<OneD, NekDouble> m_normalDiagonal;
    Array<OneD,NekDouble> m_constantProjection,m_constantComplement;
    double m_constraintMeasure = 0;
    int m_normalDirectLimit = 1024, m_normalDirectSize = 0;
    std::vector<int> m_normalDirectMap;
    std::vector<double> m_normalCholesky;
    std::vector<Array<OneD, NekDouble>> m_sampleCoordinates;
    std::vector<SolverUtils::RANSWallFrame> m_frames;
    Array<OneD, NekDouble> m_x, m_span, m_roughness, m_logRoughness, m_mask;
    Array<OneD, NekDouble> m_multiplier, m_normalTarget, m_normalError, m_projectedError;
    Array<OneD,NekDouble> m_streamDivergence,m_spanDivergence;
    Fields m_diagonal, m_rhs;
    Array<OneD, NekDouble> m_velocityMass;
    double m_beta = 0;
    mutable int m_robinNext = -1;
    mutable Fields m_robinVelocity;
    mutable Array<OneD, NekDouble> m_robinNormal, m_robinProjected;
    struct Point
    {
        SolverUtils::AWFIndependentSamples sample{};
        SolverUtils::EquilibriumAWFValues selected{};
        double kP = 0, p = 0, d = 0, epsH = 0, wallRe = 0, cutoff = 0;
        double massResidual = 0, tangentSource = 0, kSource = 0;
        double resolvedShear = 0, tangentFluxResidual = 0, kFluxResidual = 0;
    };
    std::vector<Point> m_points;
    std::ofstream m_diagnostics, m_pointDiagnostics;
};
}
#endif
