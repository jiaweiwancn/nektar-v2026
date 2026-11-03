///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SOLVERUTILS_RANS_WALL_INNER_STRESS_H
#define NEKTAR_SOLVERUTILS_RANS_WALL_INNER_STRESS_H

#include <SolverUtils/Turbulence/RANSModels.h>
#include <string>

namespace Nektar::SolverUtils
{
enum class RANSWallInnerStressModel { Linear, SaturatedR2 };
inline constexpr double RANSWallInnerSaturationRho = .8;

SOLVER_UTILS_EXPORT const char *RANSWallInnerStressName(RANSWallInnerStressModel model);
SOLVER_UTILS_EXPORT RANSWallInnerStressModel ParseRANSWallInnerStress(const std::string &name);

/** Full physical symmetric strain, turbulent mean stress F, and Reynolds
 * covariance R=(2k/3)I-F.  Molecular viscosity is deliberately not part of F.
 * The scalar auxiliary viscosity remains available to the unchanged scalar
 * diffusion and reference/pre-profile laws; it is not replaced by t/s.
 */
struct RANSWallInnerStressValues
{
    RANSWallInnerStressModel model=RANSWallInnerStressModel::Linear;
    RANSTensor strain{},turbulentStress{},reynoldsStress{};
    double k=0,auxiliaryNu=0,strainNorm=0,c=0,a=0,t=0;
    double ta=0,tc=0,radialTangent=0,tangentialTangent=0,production=0;
    bool active=false,linearResponse=true,exponentialUnderflow=false;
};

struct RANSWallInnerStressAction
{
    RANSTensor strain{},turbulentStress{},reynoldsStress{};
    double production=0;
    // At the joint k/strain or k/nuA origins this is an admissible one-sided
    // directional action, not a linear Frechet Jacobian.
    bool directionalOnly=false;
};

/** Linear is the original law, including its covariance failure at zero local
 * k with nonzero nuA*S. SaturatedR2 has fixed rho=.8, t=a below rho*c and
 * t=c-(1-rho)c exp(-(a-rho*c)/((1-rho)c)) above. At k=0 it defines F=R=0.
 * k and nuA must be nonnegative and finite; no negative extension, floor,
 * eigenvalue clipping or post-step correction is introduced.
 */
SOLVER_UTILS_EXPORT RANSWallInnerStressValues EvaluateRANSWallInnerStress(
    double k,double auxiliaryNu,const RANSTensor &gradient,
    RANSWallInnerStressModel model=RANSWallInnerStressModel::Linear);

/** Complete value action, including dnuA supplied from the physical sampled-K
 * profile and geometry chain rule. Local dk and sampled-K dependence are
 * distinct. At exact k=0/nuA=0, outward directions are rejected; joint-origin
 * limits are evaluated as actual homogeneous one-sided limits. Production
 * differentiates -R:G, including the -(2k/3) div U term.
 */
SOLVER_UTILS_EXPORT RANSWallInnerStressAction RANSWallInnerStressDirectional(
    double k,double auxiliaryNu,const RANSTensor &gradient,
    double directionK,double directionAuxiliaryNu,const RANSTensor &directionGradient,
    RANSWallInnerStressModel model=RANSWallInnerStressModel::Linear);

/** Mean stress in the reduced pressure pi=p/rho_f+2k/3.  The Linear and
 * inactive branches retain the original (nu+nuA)*S arithmetic. Active R2 uses
 * nu*S+F. There is no outer molecular term in this inner-only helper.
 */
SOLVER_UTILS_EXPORT RANSTensor RANSWallInnerMeanStress(
    const RANSWallInnerStressValues &values,double molecularNu,double reducedPressure);
} // namespace Nektar::SolverUtils
#endif
