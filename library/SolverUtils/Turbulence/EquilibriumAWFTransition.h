///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SOLVERUTILS_EQUILIBRIUM_AWF_TRANSITION_H
#define NEKTAR_SOLVERUTILS_EQUILIBRIUM_AWF_TRANSITION_H

#include <SolverUtils/Turbulence/EquilibriumAWFProfiles.h>
#include <SolverUtils/Turbulence/RANSWallNormalMap.h>

namespace Nektar::SolverUtils
{
/** Independent dimensionless indicators from the revised note. Vector entries
 * retain all Cartesian components; the caller supplies geometric tangential
 * pressure/body and convection forces separately. Nothing is divided by a
 * velocity direction or by a gate. */
struct AWFTransitionSensors
{
    RANSVector velocity{}, velocityDifference{}, pressure{}, convection{};
    double normal = 0, k = 0, epsilon = 0;
};

struct AWFTransitionInput
{
    double coordinate = 0, budgetMargin = 0, reynolds = 0;
    AWFTransitionSensors sensors;
};

struct SOLVER_UTILS_EXPORT AWFTransitionControls
{
    double start = -2, end = -1;
    AWFProfileSwitch state;
    void Validate() const;
};

/** Optional parameter direction, with all entries zero by default. Thresholds
 * are physical prescribed controls, not functions of the current flow. */
struct AWFTransitionControlDirection
{
    double start = 0, end = 0, on = 0, off = 0;
    double marginOff = 0, marginOn = 0, reynoldsOff = 0, reynoldsOn = 0;
    std::array<double, 7> weights{};
};

struct AWFTransitionValues
{
    double chi = 0, spatial = 0, state = 0, admissibility = 0;
    double marginGate = 0, reynoldsGate = 0, d2 = 0;
};

struct AWFTransitionResult
{
    AWFTransitionValues value, derivative, secondDerivative;
};

/** Analytic first and second directional derivatives along input+t*direction
 * and controls+t*controlDirection. Constant branches include both endpoints;
 * their first two derivatives are exactly zero. The helper never evaluates a
 * corrected profile/cutoff, including on an inadmissible inactive branch. */
SOLVER_UTILS_EXPORT AWFTransitionResult EvaluateEquilibriumAWFTransition(
    const AWFTransitionInput &input, const AWFTransitionControls &controls,
    const AWFTransitionInput &direction = {},
    const AWFTransitionControlDirection &controlDirection = {});

/** Physical samples all use identical current/reference evaluation maps.
 * sampledK is the inner r_K sample; kA and epsilonA are independent outer
 * samples at r_a>H. Matching samples are at H. pressure and convection are
 * the separately evaluated geometric tangent F_p and F_a, including their
 * reference values. Normal variations are derivatives of the algebraic
 * extension P=I-n*n^T; a geometric first variation satisfies n.dot(dn)=0. */
struct AWFTransitionPhysicalInput
{
    RANSVector normal{}, matchingVelocity{}, matchingReference{};
    RANSVector velocityA{}, velocityAReference{}, velocityB{}, velocityBReference{};
    RANSVector pressure{}, pressureReference{}, convection{}, convectionReference{};
    double sampledK = 0, height = 0, molecularNu = 0, Q = 0;
    double kA = 0, epsilonA = 0, epsilonAReference = 0;
};

struct AWFTransitionSensorResult
{
    AWFTransitionSensors value, derivative, secondDerivative;
};

SOLVER_UTILS_EXPORT AWFTransitionSensorResult EvaluateEquilibriumAWFTransitionSensors(
    const AWFTransitionPhysicalInput &input,
    const AWFTransitionPhysicalInput &direction = {});

struct AWFTransitionCoordinateResult
{
    double coordinate = 0, derivative = 0, secondDerivative = 0;
    RANSVector chartGradient{};
};

/** Prescribed signed wall-curve arc length from q1=0, increasing with q1.
 * Flat/plane walls allow a general constant 3D normal: at fixed q2 this is
 * q1*sqrt(1+slopeX^2). The cosine wall is extruded. Only the ground coordinates
 * are read; rho is ignored and both normal derivatives are identically zero.
 * Other geometries have no declared flat upstream reference and are rejected.
 */
SOLVER_UTILS_EXPORT double EquilibriumAWFWallCoordinate(
    const RANSWallNormalMap &wall, const RANSVector &groundChart);
SOLVER_UTILS_EXPORT AWFTransitionCoordinateResult EquilibriumAWFWallCoordinateDifferential(
    const RANSWallNormalMap &wall, const RANSVector &groundChart,
    const RANSVector &direction = {});

/** For a cosine hill, the spatial transition must finish strictly before its
 * incoming (left) foot. Flat/plane supports need no foot buffer. This verifies
 * geometry/policy, not the additional native trace-projection support buffer. */
SOLVER_UTILS_EXPORT void ValidateEquilibriumAWFTransitionSupport(
    const RANSWallNormalMap &wall, const AWFTransitionControls &controls);

/** A prescribed affine physical coordinate on a flat/plane wall. Direction
 * must be a nonzero tangent covector. No instantaneous flow enters this map. */
SOLVER_UTILS_EXPORT AWFTransitionCoordinateResult EquilibriumAWFAffineWallCoordinate(
    const RANSWallNormalMap &wall, const RANSVector &groundChart,
    const RANSVector &origin, const RANSVector &covector,
    const RANSVector &direction = {});

struct AWFTransitionProfileInput
{
    double originalDissipation = 0, correctedDissipation = 0;
    double referenceProduction = 0, production = 0;
};
struct AWFTransitionProfileValues
{
    double dissipation = 0, redistribution = 0, kSource = 0;
};
struct AWFTransitionProfileResult
{
    AWFTransitionProfileValues value, derivative, secondDerivative;
};

/** Product-rule composition with an existing profile. Profile derivatives
 * must include the current-K dependence of its implicit metric cutoff. An
 * identically inactive chi branch directly returns the baseline, without
 * inspecting correctedDissipation or referenceProduction (they may be absent).
 * This assembles the k source once, not source plus redistribution twice. */
SOLVER_UTILS_EXPORT AWFTransitionProfileResult BlendEquilibriumAWFTransition(
    const AWFTransitionResult &transition, const AWFTransitionProfileInput &profile,
    const AWFTransitionProfileInput &direction = {},
    const AWFTransitionProfileInput &secondDirection = {});
} // namespace Nektar::SolverUtils
#endif
