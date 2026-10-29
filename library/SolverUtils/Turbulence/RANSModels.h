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
enum class RANSModel { StandardKEpsilon, CLSHighRe };

struct SOLVER_UTILS_EXPORT RANSParameters
{
    double cMu = 0.09, cEpsilon1 = 1.44, cEpsilon2 = 1.92;
    double sigmaK = 1.0, sigmaEpsilon = 1.3;
    double yapCoefficient = 0.83, yapLength = 2.5;
    void Validate() const;
};

struct RANSClosureValues
{
    double cMu, eddyViscosity, strainInvariant, rotationInvariant, divergence;
    RANSTensor reynoldsStress, quadraticStress, cubicStress;
    double production, kSource, epsilonProduction, epsilonDestruction, yap;
    double epsilonSource, momentumDiffusion, kDiffusion, epsilonDiffusion;
};

SOLVER_UTILS_EXPORT double CLSCmu(double invariant);
SOLVER_UTILS_EXPORT double CLSYapFunction(double lengthRatio,
                                        const RANSParameters &parameters = {});
SOLVER_UTILS_EXPORT std::array<double, 7> CLSCubicCoefficients(double cMu);

/** grad[i][j] = dU_i/dx_j in a fixed nonrotating Cartesian frame. S and Omega
 * have no factor 1/2. The model assumes incompressibility; the raw gradient is
 * used and its divergence is reported, not silently projected away. For CLS,
 * epsilon is the transported modified dissipation in the high-Re outer limit,
 * equal to the k sink here. No damping, floors, production/realizability clips,
 * low-Re E term, or dissipation split is enabled by this evaluator.
 * resolvedNu affects ONLY diffusion. Physical AWF viscosity is separate.
 */
SOLVER_UTILS_EXPORT RANSClosureValues EvaluateRANSClosure(
    RANSModel model, double k, double epsilon, const RANSTensor &gradient,
    double wallDistance, const RANSParameters &parameters = {},
    double resolvedNu = 0.0);

// Explicit conversion for a caller choosing the finite-molecular split. It
// does not turn the high-Re evaluator into the full low-Re CLS model.
SOLVER_UTILS_EXPORT double RANSTotalDissipation(
    double k, double modifiedEpsilon, double molecularNu,
    const RANSVector &gradientK);

struct NeutralABLValues
{
    double velocity, k, epsilon, eddyViscosity, velocityGradient;
    double epsilonGradient, epsilonGradientFlux, shearStress;
};

class SOLVER_UTILS_EXPORT RANSNeutralReference
{
public:
    static RANSNeutralReference Standard(const RANSParameters &parameters = {});
    static RANSNeutralReference CLS(const RANSParameters &parameters = {});
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
};
} // namespace Nektar::SolverUtils
#endif
