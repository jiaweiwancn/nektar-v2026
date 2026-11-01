///////////////////////////////////////////////////////////////////////////////
// File: RANSVelocityCorrectionScheme.h
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_RANS_VELOCITY_CORRECTION_SCHEME_H
#define NEKTAR_RANS_VELOCITY_CORRECTION_SCHEME_H

#include <IncNavierStokesSolver/EquationSystems/VelocityCorrectionSchemeWeakPressure.h>
#include <SolverUtils/Turbulence/RANSModels.h>
#include <IncNavierStokesSolver/EquationSystems/RANSAnalyticalWallStrip.h>
#include <IncNavierStokesSolver/EquationSystems/RANSTerrainWall.h>
#include <IncNavierStokesSolver/EquationSystems/RANSMatrixFree.h>
#include <IncNavierStokesSolver/EquationSystems/RANSHistory.h>
#include <fstream>

namespace Nektar
{
/** Variable-diffusion standard/CLS k-epsilon transport. Each IMEX implicit stage
 * solves the nonlinear stress/source equations and incompressibility together.
 * Pressure is pi = p/rho + 2k/3. A pressure-defect iteration replaces the
 * constant-viscosity curl-curl extrapolation used by the legacy VCS path.
 */
class RANSVelocityCorrectionScheme : public VCSWeakPressure
{
public:
    friend class MemoryManager<RANSVelocityCorrectionScheme>;
    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &session,
        const SpatialDomains::MeshGraphSharedPtr &graph)
    {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<RANSVelocityCorrectionScheme>::AllocateSharedPtr(
                session, graph);
        p->InitObject();
        return p;
    }
    static std::string className;

protected:
    RANSVelocityCorrectionScheme(
        const LibUtilities::SessionReaderSharedPtr &session,
        const SpatialDomains::MeshGraphSharedPtr &graph);
    void v_InitObject(bool declareField = true) override;
    void v_DoInitialise(bool dumpInitialConditions = true) override;
    void v_SetInitialConditions(NekDouble time, bool dump, int domain) override;
    void v_GenerateSummary(SolverUtils::SummaryList &summary) override;
    void v_EvaluateAdvection_SetPressureBCs(
        const Array<OneD, const Array<OneD, NekDouble>> &input,
        Array<OneD, Array<OneD, NekDouble>> &output, NekDouble time) override;
    void v_SolveUnsteadyStokesSystem(
        const Array<OneD, const Array<OneD, NekDouble>> &input,
        Array<OneD, Array<OneD, NekDouble>> &output, NekDouble time,
        NekDouble lambda) override;
    void v_GetPressure(
        const Array<OneD, const Array<OneD, NekDouble>> &fields,
        Array<OneD, NekDouble> &pressure) override;
    void v_ExtraFldOutput(std::vector<Array<OneD, NekDouble>> &coefficients,
                          std::vector<std::string> &variables) override;
    bool v_PostIntegrate(int step) override;

private:
    using Fields = Array<OneD, Array<OneD, NekDouble>>;
    using Gradients = Array<OneD, Fields>;
    bool CheckState(const Fields &state, bool requirePositive = true) const;
    void UpdateLogWall(const Fields &state);
    void UpdateOutflow(const Array<OneD, NekDouble> &pressure, NekDouble time);
    void MonitorReference();
    std::string RestartSignature() const;
    void UpdateWall(const Fields &state, const Fields &predictor,
                    const Gradients &gradient,
                    const Array<OneD, NekDouble> &pressure,
                    NekDouble time, NekDouble alpha);
    void AuditStage(const Fields &state, const Fields &predictor,
                    const Fields &diffusion, const Fields &gain, const Fields &loss,
                    const Array<OneD, NekDouble> &pressure,
                    NekDouble time, NekDouble alpha, int phase);
    void Closure(const Fields &state, Gradients &gradient, Fields &diffusion,
                 Fields &gain, Fields &loss, NekDouble time);
    void SolveField(int field, const Array<OneD, NekDouble> &rhs,
                    const Array<OneD, NekDouble> &diffusion,
                    const Array<OneD, NekDouble> &mass,
                    Array<OneD, NekDouble> &output, NekDouble alpha,
                    const Array<OneD, NekDouble> *newtonTrial = nullptr);
    NekDouble CorrectPressure(const Fields &velocity,
                              const Array<OneD, NekDouble> &mu,
                              Array<OneD, NekDouble> &pressure,
                              NekDouble alpha, NekDouble &divergence,
                              NekDouble &weakDivergence);
    NekDouble Maximum(NekDouble value) const;
    void RequireEverywhere(bool valid, const std::string &message) const;

    SolverUtils::RANSParameters m_parameters;
    SolverUtils::RANSModel m_model;
    Array<OneD, NekDouble> m_wallDistance;
    Gradients m_nonlinearStress;
    std::shared_ptr<RANSWallBoundary> m_wallStrip;
    SolverUtils::RANSWallGeometry m_wallGeometry;
    bool m_terrainWall = false;
    std::shared_ptr<RANSMatrixFree> m_matrixFree;
    std::shared_ptr<RANSIMEX> m_historyScheme;
    size_t m_historyStepOffset = 0;
    Fields m_importedCoefficients;
    int m_linearIterations = 0;
    NekDouble m_wallOrigin;
    std::shared_ptr<const SolverUtils::RANSNeutralReference> m_reference;
    MultiRegions::ExpListSharedPtr m_pressureMass;
    MultiRegions::ExpListSharedPtr m_pressurePoisson;
    NekDouble m_pressureInverseAlpha = 0;
    int m_k, m_epsilon, m_maxIterations;
    bool m_useNewton = false;
    bool m_preconditionedNewton = false;
    bool m_pointwiseNewtonScalars = false;
    NekDouble m_newtonDifference = 1e-6;
    int m_newtonKrylovDimension = 240;
    bool m_nonlinearTrial = false, m_freezeCmuBranch = false;
    std::vector<int> m_cmuBranches;
    NekDouble m_resolvedNu, m_tolerance, m_pressureTolerance, m_pressureRelaxation, m_relaxation;
    NekDouble m_minimumK, m_minimumEpsilon, m_z0, m_heightShift;
    NekDouble m_gradDiv = 0;
    bool m_implicitAdvection = false;
    // Optional discrete continuation for an experimental flat precursor.
    // It adds a fading initial weak residual; no field is reset or clipped.
    NekDouble m_startupDuration = 0;
    Fields m_startupLoad;
    NekDouble StartupWeight(NekDouble time) const;
    bool m_auditStage = false;
    bool m_monitorReference = false, m_hasOutflow = false;
    std::vector<int> m_outflows;
    std::string m_wallModel;
    std::ofstream m_diagnostics;
    std::ofstream m_residualDiagnostics;
    std::ofstream m_referenceDiagnostics;
    static std::string solverTypeLookupId;
};
} // namespace Nektar
#endif
