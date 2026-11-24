// SPDX-License-Identifier: MIT
#ifndef NEKTAR_WJ_NATIVE_FLOW_TANGENT_H
#define NEKTAR_WJ_NATIVE_FLOW_TANGENT_H
#include <IncNavierStokesSolver/EquationSystems/WJNativeRules.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeWallV2.h>
#include <array>

namespace Nektar::WJNative
{
/** Exact partial derivatives of the approved V tensor at fixed k, epsilon,
 * auxiliary viscosity, radius, H and physical nu. Input columns are ordered
 * G_xx,G_xz,G_zx,G_zz; each output remains a full 3D tensor. The WJ strain still
 * subtracts trace(G)/3, including its nonzero yy component under weak divergence.
 */
struct FrozenPlanarVelocityTangent
{
    std::array<Tensor,4> dVdG{};
};
struct OwnedPlanarVelocityTangents
{
    int geometry=-1;
    std::vector<FrozenPlanarVelocityTangent> points;
};

FrozenPlanarVelocityTangent EvaluatePlanarVelocityTangent(double k,double epsilon,
    const Tensor&gradient,double nu,double auxiliary,double r,double H,
    ConstitutiveBackend backend=ConstitutiveBackend::Historical,
    WallEvaluationStats*statistics=nullptr, const Vector&normal={0,0,1});

/** Full signed canonical [u,v] Jacobian of int grad(test):V, size 2Nu by 2Nu.
 * Samples must match each supplied owned rule point in its original order.
 * NativeRule B/Dx/Dz already contain signed assembly factors; no second sign
 * multiplication is performed. Essential reaction rows and all cross blocks
 * are retained. This is generally nonsymmetric/indefinite: no positivity
 * modification, symmetrization, pressure, advection, mass, boundary, or scalar
 * derivative is added. The caller owns scalar/Gamma coupling and admissibility.
 */
Sparse AssemblePlanarVelocityTangent(const Rules&rules,int velocitySize,
    const std::vector<OwnedPlanarVelocityTangents>&samples);
}
#endif
