///////////////////////////////////////////////////////////////////////////////
// File: RANSModels.h
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SOLVERUTILS_RANS_MODELS_H
#define NEKTAR_SOLVERUTILS_RANS_MODELS_H

#include <SolverUtils/SolverUtilsDeclspec.h>
#include <array>
#include <string>

namespace Nektar::SolverUtils
{
using RANSVector = std::array<double, 3>;
using RANSTensor = std::array<RANSVector, 3>;
enum class RANSModel { StandardKEpsilon, CLSHighRe, A93, A95, E94E95 };
enum class RANSKappaPolicy { Natural, Prescribed };
enum class RANSWallFormulation { EPAWF, AcceptedBudgetAWF };
SOLVER_UTILS_EXPORT bool IsNewRANSModel(RANSModel model);
SOLVER_UTILS_EXPORT const char *RANSModelName(RANSModel model);
SOLVER_UTILS_EXPORT RANSModel ParseRANSModel(const std::string &name);

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
 * resolvedNu affects ONLY diffusion for legacy models. For A93/A95/E94E95,
 * resolve parameters using RANSResolvedModel: tensors are physical half-scaled,
 * invariants/nonlinear tensors use deviatoric strain, production contracts full
 * R with raw A, and resolvedNu must be zero (high-Re outer equations).
 */
SOLVER_UTILS_EXPORT RANSClosureValues EvaluateRANSClosure(
    RANSModel model, double k, double epsilon, const RANSTensor &gradient,
    double wallDistance, const RANSParameters &parameters = {},
    double resolvedNu = 0.0);

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
};

/** Immutable resolved physics descriptor. The historical type name is retained
 * for CLS source/ABI naming; RANSResolvedModel is the preferred new-model name.
 * New models use total epsilon, high-Re outer diffusion, finite physical inner
 * viscosity, DeviatoricStrainV1 and no constitutive/source regularization.
 */
class SOLVER_UTILS_EXPORT RANSNeutralReference
{
public:
    static RANSNeutralReference Resolve(RANSModel model, RANSKappaPolicy policy,
        double prescribedKappa, RANSWallFormulation wall, double molecularNu);
    bool NewModel() const { return IsNewRANSModel(m_model); }
    RANSKappaPolicy KappaPolicy() const { return m_policy; }
    RANSWallFormulation WallFormulation() const { return m_wall; }
    double MolecularNu() const { return m_molecularNu; }
    std::string Signature() const;
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
    RANSKappaPolicy m_policy = RANSKappaPolicy::Natural;
    RANSWallFormulation m_wall = RANSWallFormulation::EPAWF;
    double m_molecularNu = 0;

};
using RANSResolvedModel = RANSNeutralReference;
SOLVER_UTILS_EXPORT RANSClosureValues EvaluateRANSClosure(
    const RANSResolvedModel &model, double k, double epsilon,
    const RANSTensor &gradient, double wallDistance);
} // namespace Nektar::SolverUtils
#endif
