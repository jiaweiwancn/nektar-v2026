///////////////////////////////////////////////////////////////////////////////
// File: RANSModels.h
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SOLVERUTILS_RANS_MODELS_H
#define NEKTAR_SOLVERUTILS_RANS_MODELS_H

#include <SolverUtils/SolverUtilsDeclspec.h>
#include <array>

namespace Nektar::SolverUtils
{
using RANSVector = std::array<double, 3>;
using RANSTensor = std::array<RANSVector, 3>;
enum class RANSModel { StandardKEpsilon, CLSHighRe, WJHighRe };
struct WJParameters;

struct SOLVER_UTILS_EXPORT RANSParameters
{
    double cMu = 0.09, cEpsilon1 = 1.44, cEpsilon2 = 1.92;
    double sigmaK = 1.0, sigmaEpsilon = 1.3;
    double yapCoefficient = 0.83, yapLength = 2.5;
    void Validate() const;
};

struct WJDiagnostics
{
    // q is the tensor coefficient denominator Q, not the neutral shear qe.
    double n = 0, s = 0, w = 0, i4 = 0, i5 = 0, q = 0;
    double rootResidual = 0, rootDerivative = 0;
    std::array<double, 5> beta{}; // beta_1, beta_3, beta_4, beta_6, beta_9
};

struct RANSClosureValues
{
    double cMu, eddyViscosity, strainInvariant, rotationInvariant, divergence;
    RANSTensor reynoldsStress, quadraticStress, cubicStress;
    double production, kSource, epsilonProduction, epsilonDestruction, yap;
    double epsilonSource, momentumDiffusion, kDiffusion, epsilonDiffusion;
    // Generic stress split: R = 2k I/3 + deviatoricStress,
    // deviatoricStress = -2 nu_eff dev(sym(G)) + stressRemainder for WJ.
    // quadraticStress/cubicStress are CLS-specific and remain zero for WJ.
    RANSTensor anisotropy{}, deviatoricStress{}, stressRemainder{};
    double scalarViscosity = 0, productionDev = 0;
    WJDiagnostics wj{};
};

SOLVER_UTILS_EXPORT double CLSCmu(double invariant);
SOLVER_UTILS_EXPORT double CLSYapFunction(double lengthRatio,
                                        const RANSParameters &parameters = {});
SOLVER_UTILS_EXPORT std::array<double, 7> CLSCubicCoefficients(double cMu);

/** Legacy Standard/CLS API: grad[i][j] = dU_i/dx_j in a fixed nonrotating Cartesian frame. S and Omega
 * have no factor 1/2. The model assumes incompressibility; the raw gradient is
 * used and its divergence is reported, not silently projected away. For CLS,
 * epsilon is the transported modified dissipation in the high-Re outer limit,
 * equal to the k sink here. No damping, floors, production/realizability clips,
 * low-Re E term, or dissipation split is enabled by this evaluator.
 * resolvedNu affects ONLY diffusion. Physical AWF viscosity is separate.
 */
SOLVER_UTILS_EXPORT RANSClosureValues EvaluateRANSClosure(
    RANSModel model, double k, double epsilon, const RANSTensor &gradient,
    double wallDistance, const RANSParameters &parameters,
    double resolvedNu = 0.0);

// Model-aware omitted-parameters dispatch. WJ uses its distinct defaults.
SOLVER_UTILS_EXPORT RANSClosureValues EvaluateRANSClosure(
    RANSModel model, double k, double epsilon, const RANSTensor &gradient,
    double wallDistance);

/** Generalized-Jacobian probe only: choose the active max(strain,rotation)
 * branch at the base iterate, while differentiating that branch's invariant.
 * The physical residual and all reported stresses use EvaluateRANSClosure.
 */
SOLVER_UTILS_EXPORT RANSClosureValues EvaluateRANSClosureLinearization(
    RANSModel model, double k, double epsilon, const RANSTensor &gradient,
    double wallDistance, const RANSParameters &parameters, double resolvedNu,
    bool strainBranch);

// Explicit conversion for a caller choosing the finite-molecular split. It
// does not turn the high-Re evaluator into the full low-Re CLS model.
SOLVER_UTILS_EXPORT double RANSTotalDissipation(
    double k, double modifiedEpsilon, double molecularNu,
    const RANSVector &gradientK);

struct NeutralABLValues
{
    double velocity, k, epsilon, eddyViscosity, velocityGradient;
    double epsilonGradient, epsilonGradientFlux, shearStress;
    RANSTensor covariance{}; // Cartesian (x,y,z), retaining transverse yy
};

class SOLVER_UTILS_EXPORT RANSNeutralReference
{
public:
    static RANSNeutralReference Standard(const RANSParameters &parameters = {});
    static RANSNeutralReference CLS(const RANSParameters &parameters = {});
    static RANSNeutralReference WJ();
    static RANSNeutralReference WJ(const WJParameters &parameters);
    RANSModel Model() const { return m_model; }
    double Q() const { return m_q; }
    double Kappa() const { return m_kappa; }
    double Cmu() const { return m_cmu; }
    const RANSParameters &Parameters() const { return m_parameters; }

    // Unshifted geometric height; anchored to an actual matching value.
    NeutralABLValues AtHeight(double height, double frictionVelocity,
                              double matchingHeight, double matchingVelocity) const;
    // Richards--Norris shifted-origin option; prohibited for CLS/Yap geometry.
    NeutralABLValues StandardShifted(double physicalHeight,
                                     double frictionVelocity, double z0) const;
    std::array<double, 2> StandardLogWallTraction(
        double kP, const std::array<double, 2> &velocityP,
        double sampleHeight, double z0, double heightShift = 0.0) const;

private:
    RANSNeutralReference(RANSModel model, const RANSParameters &parameters,
                         double q, double kappa, double cmu);
    RANSModel m_model;
    RANSParameters m_parameters;
    double m_q, m_kappa, m_cmu;
    double m_scalarRatio = 1;
};
} // namespace Nektar::SolverUtils
#endif
