///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SOLVERUTILS_WJ_EARSM_H
#define NEKTAR_SOLVERUTILS_WJ_EARSM_H

#include <SolverUtils/Turbulence/RANSModels.h>

namespace Nektar::SolverUtils
{
/** Separate high-Re WJ completion: epsilon is transported directly. Cd is a
 * state-independent scalar coefficient, never the effective stress Cmu.
 * Only the selected c1/c2 tensor law is supported. No legacy knobs are accepted.
 */
struct SOLVER_UTILS_EXPORT WJParameters
{
    static constexpr double Cstar = 572.0/6561.0;
    double c1 = 1.8, c2 = 5.0/9.0, cd = Cstar;
    double cEpsilon1 = 1.44, cEpsilon2 = 1.82;
    double sigmaK = 1, sigmaEpsilon = 1.3;
    static WJParameters Defaults() { return {}; }
    void Validate() const;
};

enum class WJGradientPolicy
{
    Incompressible,
    // Off-constraint completion: dev strain for covariance, raw stress work
    // for scalar sources. Agrees with the notes on the incompressible manifold.
    DevStrainRawWorkV1
};

SOLVER_UTILS_EXPORT RANSClosureValues EvaluateWJClosure(double k, double epsilon,
    const RANSTensor &rawGradient, const WJParameters &parameters = WJParameters::Defaults(),
    WJGradientPolicy policy = WJGradientPolicy::Incompressible, double outerNu = 0);

/** Common exact directional action with respect to k, epsilon and raw G.
 * All entries, including diagnostics and coefficients, are directional values.
 * N uses the implicit physical-root derivative; no iteration is differentiated.
 */
SOLVER_UTILS_EXPORT RANSClosureValues ApplyWJClosureDirection(double k, double epsilon,
    const RANSTensor &rawGradient, double dk, double dEpsilon,
    const RANSTensor &dRawGradient,
    const WJParameters &parameters = WJParameters::Defaults(),
    WJGradientPolicy policy = WJGradientPolicy::Incompressible, double outerNu = 0);

/** Kernel-level factory/adapter only. Native PDE selection remains disabled. */
class SOLVER_UTILS_EXPORT RANSClosureKernel
{
public:
    static RANSClosureKernel WJ(
        const WJParameters &parameters = WJParameters::Defaults(),
        WJGradientPolicy policy = WJGradientPolicy::Incompressible);
    RANSClosureValues Evaluate(double k, double epsilon, const RANSTensor &rawG,
                              double outerNu = 0) const;
    RANSClosureValues Direction(double k, double epsilon, const RANSTensor &rawG,
        double dk, double dEpsilon, const RANSTensor &dRawG, double outerNu = 0) const;
private:
    RANSClosureKernel(const WJParameters &parameters, WJGradientPolicy policy);
    WJParameters m_parameters;
    WJGradientPolicy m_policy;
};
}
#endif
