///////////////////////////////////////////////////////////////////////////////
// File: EquilibriumAWFProfiles.h
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SOLVERUTILS_EQUILIBRIUM_AWF_PROFILES_H
#define NEKTAR_SOLVERUTILS_EQUILIBRIUM_AWF_PROFILES_H

#include <SolverUtils/Turbulence/AnalyticalWallFunction.h>
#include <functional>
#include <vector>

namespace Nektar::SolverUtils
{
/** Physical normal distances. sampleHeight samples k only; tailStart is an
 * independent dissipation-tail requirement. Neither is an FV P/N location. */
struct SOLVER_UTILS_EXPORT AWFProfileGeometry
{
    double molecularNu, height, roughness, sampleHeight, tailStart;
    void Validate() const;
    double SampleK(const std::function<double(double)> &field) const;
};

/** B=grad_Gamma(n), n into the fluid: J=(1+kappa1*r)(1+kappa2*r).
 * This validates one ray's local metric, not global normal-tube injectivity. */
struct SOLVER_UTILS_EXPORT AWFNormalMetric
{
    double kappa1 = 0, kappa2 = 0, minimumStretch = .1;
    void Validate(double height) const;
    double Jacobian(double r) const;
};

struct SOLVER_UTILS_EXPORT AWFProfileSwitch
{
    double on = .02, off = .1, marginOff = .05, marginOn = .1;
    double reynoldsOff = 100, reynoldsOn = 200;
    std::array<double, 7> weights{1, 1, 1, 1, 1, 1, 1};
    void Validate() const;
};

/** Dimensionless independent indicators in order U, DeltaU, normal, k,
 * epsilon, pressure, convection. Physical sampling/projection is the caller's
 * responsibility. This scalar weight is extended constantly along a ray. */
SOLVER_UTILS_EXPORT double EquilibriumAWFProfileWeight(
    double spatialWeight, const std::array<double, 7> &sensors,
    double budgetMargin, double reynolds, const AWFProfileSwitch &controls = {});

/** Wall-source formulation, independent of the outer constitutive variant. */
enum class AWFProfileFormulation { LegacyEPAWF, CLSCalibratedBudgetAWF, AcceptedBudgetAWF };

/** Separate auxiliary-wall parameters derived from the selected CLS reference.
 * Keeping this explicit prevents factories from silently retaining the legacy
 * 0.09/2.55 wall coefficients for the calibrated formulation. */
SOLVER_UTILS_EXPORT AWFParameters CalibratedAWFParameters(const RANSNeutralReference &model);

enum class AWFCurvedCapability { FlatQualifiedV1, AcceptedVerificationV1 };

// Physical sink policy, independent of closure and numerical geometry version.
// The displaced policy retains the original physical wall and resistance.
enum class AWFDisplacementPolicy { UnshiftedV1, SamePhysicalWallDisplacedAcceptedV1 };

struct AWFProfileOptions
{
    AWFCurvedCapability curvedCapability = AWFCurvedCapability::FlatQualifiedV1;
    AWFDisplacementPolicy displacementPolicy = AWFDisplacementPolicy::UnshiftedV1;
    double aerodynamicDisplacement = 0;
    double requestedWeight = 1, minimumReynolds = 100;
    int quadratureOrder = 16;
    AWFFailurePolicy failurePolicy = AWFFailurePolicy::RequirePreservation;
    // Additional physical mesh interfaces: the cutoff and every source load
    // must integrate over these same panels when the wall spans elements.
    std::vector<double> physicalBreakpoints;
    AWFProfileFormulation formulation = AWFProfileFormulation::LegacyEPAWF;
};

struct AWFProfileQuadraturePoint
{
    double r, weight, jacobian; // weight is dr; weight*jacobian is physical volume/wall area
    // Generating physical panel, retained before a Gauss coordinate is rounded.
    // Callers mapping to native elements must use the panel, not relocate r.
    double panelLower = 0, panelUpper = 0;
};

struct AWFProfileReference
{
    double frictionVelocity = 0, traction = 0, matchingVelocity = 0, numerator = 0;
    double productionBudget = 0, minimumBudget = 0, margin = 0, reynolds = 0;
    double cutoff = 0, requestedWeight = 0, weight = 0, budgetResidual = 0;
    int cutoffIterations = 0;
    AWFFallbackReason reason = AWFFallbackReason::SwitchInactive;
};

struct AWFProfileValues
{
    double production, referenceProduction, dissipation, redistribution, kSource;
};

/** Immutable local profile/load operator from the revised note.
 * Retains original AWF momentum and correlations; supplies E_EP and the signed
 * zero-budget redistribution. The simplified source is assembled once.
 * A curved/full-strain patch supplies its actual production to AtProduction;
 * flat parallel Production is only the scalar-kernel reduction.
 *
 * Source loads are dual loads per unit physical ground area, with all supplied
 * basis rows evaluated at Quadrature(). Tangential transfer, a multiregion
 * field layout, and the global patch solve are separate solver responsibilities.
 */
class SOLVER_UTILS_EXPORT EquilibriumAWFProfiles
{
public:
    EquilibriumAWFProfiles(const AWFProfileGeometry &geometry, double sampledK,
        const RANSNeutralReference &model, const AWFNormalMetric &metric = {},
        const AWFProfileOptions &options = {}, const AWFParameters &wall = {});
    const AWFProfileGeometry &Geometry() const { return m_geometry; }
    const AWFNormalMetric &Metric() const { return m_metric; }
    const AWFProfileOptions &Options() const { return m_options; }
    AWFProfileFormulation Formulation() const { return m_options.formulation; }
    bool AcceptedAWF() const
    { return Formulation() == AWFProfileFormulation::CLSCalibratedBudgetAWF ||
             Formulation() == AWFProfileFormulation::AcceptedBudgetAWF; }
    const AnalyticalWallFunction &Original() const { return m_original; }
    const AWFProfileReference &Reference() const { return m_reference; }
    const std::vector<AWFProfileQuadraturePoint> &Quadrature() const { return m_quadrature; }
    // Read-only integration panels, retaining this profile's solved cutoff.
    std::vector<AWFProfileQuadraturePoint> IntegrationQuadrature(double upper,
        const std::vector<double> &additionalBreakpoints={}) const;
    double K() const { return m_original.K(); }
    double Production(double r, const AWFVector &velocity, const AWFVector &forcing) const;
    double ReferenceProduction(double r) const;
    AWFProfileValues AtProduction(double r, double production) const;
    AWFProfileValues Evaluate(double r, const AWFVector &velocity, const AWFVector &forcing) const;
    double Dissipation(double r) const;
    double RedistributionFlux(double r) const; // includes chi, normal component in physical coordinates
    /** Stationary reference only; never a field-restoration operation.
     * For accepted AWF, ReferenceKFlux is Q_J=J D_k k'_ref, so the physical
     * upward normal diffusive flux is -Q_J/J. With J=1 it is appendix Q_0.
     * Legacy EP references retain constant k and zero ordinary diffusion.
     * A nonpositive auxiliary reference is returned without clipping: cutoff
     * applicability and qualification of that reference are separate checks. */
    double ReferenceK(double r) const;
    double ReferenceKGradient(double r) const;
    double ReferenceKFlux(double r) const;
    double ReferenceKDiffusivity(double r) const;
    // Metric-weighted corrected-dissipation integral at an independent cutoff.
    double DissipationBudget(double cutoff) const;
    // Reuse an already qualified cutoff when only the ray-constant weight
    // changes. A zero weight is constructed separately from the baseline,
    // whose quadrature must not require or retain an inadmissible cutoff.
    EquilibriumAWFProfiles Reweighted(double weight) const;
    std::vector<double> RayLoads(const std::vector<std::vector<double>> &basis,
                                 const std::vector<double> &samples) const;
    std::vector<double> SourceLoads(const std::vector<std::vector<double>> &basis,
        const AWFVector &velocity, const AWFVector &forcing) const;

private:
    void CheckHeight(double r) const;
    void GaussRule();
    std::vector<AWFProfileQuadraturePoint> MakeQuadrature(double upper, double cutoff,
                                                       double lower = 0) const;
    double Integrate(const std::function<double(double)> &function,
                     double upper, double cutoff) const;
    double Corrected(double r, double cutoff) const;
    double CorrectedIntegral(double upper, double cutoff) const;
    double BudgetDifference(double cutoff) const;
    AWFProfileGeometry m_geometry;
    AWFNormalMetric m_metric;
    AWFProfileOptions m_options;
    double m_sigmaK = 1;
    AnalyticalWallFunction m_original;
    AWFProfileReference m_reference;
    std::vector<double> m_gaussPoints, m_gaussWeights;
    std::vector<AWFProfileQuadraturePoint> m_quadrature;
};
}
#endif
