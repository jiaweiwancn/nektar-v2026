// SPDX-License-Identifier: MIT
#ifndef NEKTAR_WJ_NATIVE_WALL_V2_H
#define NEKTAR_WJ_NATIVE_WALL_V2_H
#include <IncNavierStokesSolver/EquationSystems/WJNativeFlow.h>
#include <cstdint>
namespace Nektar::WJNative
{
inline constexpr const char*WallModelV2="FiniteNuIsotropicWallFullWJCollarV2";
enum class ConstitutiveBackend { Historical, FastFull3DV1 };
struct WallEvaluationStats
{
    std::uint64_t values=0,directions=0,fallbacks=0,rootIterations=0;
    double maximumRootResidual=0;
};
WallStress AdaptAmended(double k,double epsilon,const Tensor&gradient,
    const Vector&normal,double nu,double auxiliary,double r,double H,
    ConstitutiveBackend backend=ConstitutiveBackend::Historical,WallEvaluationStats*stats=nullptr);
/** Exact direction at fixed unit normal and H, including chi's r derivative.
 * Physical nu is fixed. Auxiliary's direction includes its declared K/r law.
 */
WallStress ApplyAmendedDirection(double k,double epsilon,const Tensor&gradient,
    const Vector&normal,double nu,double auxiliary,double r,double H,
    double dk,double dEpsilon,const Tensor&dGradient,double dAuxiliary,double dr=0,
    ConstitutiveBackend backend=ConstitutiveBackend::Historical,WallEvaluationStats*stats=nullptr);
}
#endif
