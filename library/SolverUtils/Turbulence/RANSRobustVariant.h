///////////////////////////////////////////////////////////////////////////////
// File: RANSRobustVariant.h
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SOLVERUTILS_RANS_ROBUST_VARIANT_H
#define NEKTAR_SOLVERUTILS_RANS_ROBUST_VARIANT_H

#include <SolverUtils/Turbulence/RANSModels.h>
#include <vector>

namespace Nektar::SolverUtils
{
/** Explicit constitutive variant. beta scales every active cubic term by the
 * same constant. gamma has viscosity units and is prescribed independently of
 * the evolving state; referenceStrain is a symmetric, trace-free physical
 * Cartesian tensor. No invariant cap or continuation of C(q) is introduced.
 */
struct SOLVER_UTILS_EXPORT RANSRobustVariant
{
    double beta = 1, gamma = 0;
    RANSTensor referenceStrain{};
    void Validate() const;
    bool Original() const { return beta == 1 && gamma == 0; }
};

/** Calls the original kernel, then applies (beta-1)*R_cubic-gamma*(S-S_ref).
 * Production, kSource, epsilonProduction and epsilonSource use that same
 * modified stress. The original diffusion scalars and all other coefficients
 * are retained: the gamma stress must not be added a second time as diffusion.
 * beta=1,gamma=0 returns the original result without further arithmetic.
 */
SOLVER_UTILS_EXPORT RANSClosureValues EvaluateRANSRobustClosure(
    RANSModel model, double k, double epsilon, const RANSTensor &gradient,
    double wallDistance, const RANSParameters &parameters = {},
    double resolvedNu = 0, const RANSRobustVariant &variant = {});

/** Formal branch evaluator for differentiating the selected original invariant
 * branch. beta, gamma and S_ref remain prescribed. This is not an adaptive
 * coefficient or a replacement physical max(strain,rotation) residual.
 */
SOLVER_UTILS_EXPORT RANSClosureValues EvaluateRANSRobustClosureLinearization(
    RANSModel model, double k, double epsilon, const RANSTensor &gradient,
    double wallDistance, const RANSParameters &parameters, double resolvedNu,
    bool strainBranch, const RANSRobustVariant &variant = {});

struct SOLVER_UTILS_EXPORT RANSCappedDestructionControls
{
    double theta = .2, delta = .1;
    bool smooth = false; // false: selected piecewise linear limiter; true: C2 quintic
    double budgetRelativeTolerance = 1e-10;
    void Validate() const;
};

/** Positive physical metric weights are J*dr (or their positive quadrature
 * counterparts). preE>0, preD>=0 and their integrated budgets must agree to
 * the declared tolerance. Transported k>=0; sampledK>0. Negative k is rejected,
 * not extended or clipped. Corrected-cutoff construction belongs to the caller;
 * an inactive cutoff need never be evaluated to supply preE=preD=E_A.
 */
struct RANSCappedDestructionInput
{
    double sampledK = 0;
    std::vector<double> preE, preD, k, metricWeights;
};

struct RANSCappedDestructionPoint
{
    double limiter = 0, dissipation = 0, destruction = 0, redistribution = 0;
    double rate = 0, rateBound = 0;
};

struct RANSCappedDestructionValues
{
    double sampledK = 0, budgetE = 0, budgetD = 0, budgetMismatch = 0;
    double b = 0, m = 0, lambda = 0;
    double dissipationBudget = 0, destructionBudget = 0, redistributionBudget = 0;
    bool capped = false;
    std::vector<double> limiter, dissipation, destruction, redistribution, rate, rateBound;
};

/** The common budget is evaluated from preD after checking its equality with
 * preE. This preserves the all-ray limiter plateau exactly, even when the
 * existing cutoff solve leaves a roundoff-sized pre-budget mismatch. The
 * mismatch is exposed, never repaired by altering either profile. Consequently
 * the output redistribution budget equals lambda*(budgetE-budgetD) to roundoff.
 */
SOLVER_UTILS_EXPORT RANSCappedDestructionValues EvaluateRANSCappedDestruction(
    const RANSCappedDestructionInput &input,
    const RANSCappedDestructionControls &controls = {});

/** Evaluate a new location (including H) using the stored ray normalization.
 * Its local k and pre-profiles must be supplied at that location. E_pos=lambda
 * E_pre, not a local-k rescaling; the source is actual production-D_pos once.
 */
SOLVER_UTILS_EXPORT RANSCappedDestructionPoint EvaluateRANSCappedDestructionPoint(
    double k, double sampledK, double preE, double preD, double b, double m,
    const RANSCappedDestructionControls &controls = {});

SOLVER_UTILS_EXPORT RANSCappedDestructionPoint EvaluateRANSCappedDestructionPoint(
    double preE, double preD, double localK, double sampledK,
    const RANSCappedDestructionValues &normalization,
    const RANSCappedDestructionControls &controls = {});

/** Full directional input: an empty vector denotes zero in that block. The
 * controls remain fixed. All metric/profile/sampled-K dependence is included.
 */
struct RANSCappedDestructionDirection
{
    double sampledK = 0;
    std::vector<double> preE, preD, k, metricWeights;
};

/** Returns directional derivatives in the numeric fields of Values (capped is
 * the base-state flag). At limiter/cap kinks this is the actual one-sided
 * directional derivative: min/max select according to the directional sign.
 * A direction taking an exactly zero k out of the admissible cone is rejected.
 * The linear limiter rate at k=0 uses its nonnegative one-sided limit; the C2
 * limiter rate and its first k derivative both vanish there.
 */
SOLVER_UTILS_EXPORT RANSCappedDestructionValues RANSCappedDestructionDirectional(
    const RANSCappedDestructionInput &input,
    const RANSCappedDestructionDirection &direction,
    const RANSCappedDestructionControls &controls = {});
} // namespace Nektar::SolverUtils
#endif
