// SPDX-License-Identifier: MIT
#ifndef NEKTAR_WJ_NATIVE_CONSTITUTIVE_H
#define NEKTAR_WJ_NATIVE_CONSTITUTIVE_H
#include <SolverUtils/Turbulence/WJEARSM.h>

namespace Nektar::WJNative
{
enum class FastWJFallback
{
    None, ExtremeScale, InvariantRange, UnverifiedRoot, NonfiniteOutput, ViscosityUnderflow
};
struct FastWJStatus
{
    FastWJFallback fallback=FastWJFallback::None;
    int rootIterations=0;
    double normalizedRootResidual=0, relativeRootCondition=0;
    bool UsedFallback()const{return fallback!=FastWJFallback::None;}
};
struct FastWJPair
{
    SolverUtils::RANSClosureValues value{},direction{};
    FastWJStatus status;
};

/** Mathematically identical WJ-1, full 3D DevStrainRawWorkV1, outerNu=0.
 * The selected N>=9/5 root is verified before use. Extreme or unverified
 * double arithmetic delegates to the unchanged robust long-double kernel.
 * No state/input clipping, alternate root, scalar coefficient substitution,
 * or planar strain shortcut is permitted. Status is caller-owned, so the
 * caller can count every fallback without a mutable global kernel cache.
 */
SolverUtils::RANSClosureValues EvaluateFastWJClosure(double k,double epsilon,
    const SolverUtils::RANSTensor &gradient,const SolverUtils::WJParameters &parameters={},
    FastWJStatus *status=nullptr);
SolverUtils::RANSClosureValues ApplyFastWJClosureDirection(double k,double epsilon,
    const SolverUtils::RANSTensor &gradient,double dk,double dEpsilon,
    const SolverUtils::RANSTensor &dGradient,const SolverUtils::WJParameters &parameters={},
    FastWJStatus *status=nullptr);
FastWJPair EvaluateFastWJValueAndDirection(double k,double epsilon,
    const SolverUtils::RANSTensor &gradient,double dk,double dEpsilon,
    const SolverUtils::RANSTensor &dGradient,const SolverUtils::WJParameters &parameters={});
const char *FastWJFallbackName(FastWJFallback reason);
}
#endif
