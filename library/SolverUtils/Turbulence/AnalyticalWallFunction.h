///////////////////////////////////////////////////////////////////////////////
// File: AnalyticalWallFunction.h
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SOLVERUTILS_ANALYTICAL_WALL_FUNCTION_H
#define NEKTAR_SOLVERUTILS_ANALYTICAL_WALL_FUNCTION_H

#include <SolverUtils/Turbulence/RANSModels.h>
#include <array>

namespace Nektar::SolverUtils
{
using AWFVector = std::array<double, 2>;
struct SOLVER_UTILS_EXPORT AWFGeometry
{
    double molecularNu, height, roughness, pointP;
    void Validate() const;
};
struct SOLVER_UTILS_EXPORT AWFParameters
{
    double cMu = 0.09, cLength = 2.55, smoothIntercept = 10.7;
    void Validate() const;
};
struct AWFIntegrals
{
    double j0, j1, j2, productionWeight;
};
struct AWFValues
{
    AWFVector traction;
    double meanProduction, meanDissipation, dissipationP, dissipationH;
};

/** Published smooth/rough AWF, per unit density in a fixed tangential frame.
 * Positive traction is the note's shear convention; force on the fluid has the
 * opposite sign. Molecular viscosity, K, and geometry are explicit and fixed
 * for this local evaluation. Reconstruct the object when current K changes.
 */
class SOLVER_UTILS_EXPORT AnalyticalWallFunction
{
public:
    AnalyticalWallFunction(const AWFGeometry &geometry, double kP,
                           const AWFParameters &parameters = {});
    const AWFGeometry &Geometry() const { return m_geometry; }
    const AWFParameters &Parameters() const { return m_parameters; }
    double K() const { return m_k; }
    double Intercept() const;
    double ViscosityInterceptHeight() const;
    double DissipationCutoff() const;
    double EddyViscosity(double z) const;
    double TotalViscosity(double z) const;
    const AWFIntegrals &Integrals() const { return m_integrals; }
    AWFIntegrals Integrals(double z) const;
    double OriginalDissipation(double z) const;
    double MeanOriginalDissipation() const;
    AWFValues Evaluate(const AWFVector &matchingVelocity,
                       const AWFVector &forcing) const;
    AWFVector Velocity(double z, const AWFVector &matchingVelocity,
                        const AWFVector &forcing) const;
    AWFVector VelocityGradient(double z, const AWFVector &matchingVelocity,
                                const AWFVector &forcing) const;

private:
    void CheckHeight(double z) const;
    std::array<long double, 3> ProductionMoments(long double left,
                                                long double right) const;
    std::array<long double, 2> Traction(const AWFVector &velocity,
                                       const AWFVector &forcing) const;
    AWFGeometry m_geometry;
    AWFParameters m_parameters;
    double m_k;
    long double m_intercept, m_zv, m_slope, m_ze;
    AWFIntegrals m_integrals;
};

enum class AWFFallbackReason
{
    Active, SwitchInactive, CutoffInadmissible, OutsideHighReRegime,
    UnrepresentableCutoff
};
enum class AWFFailurePolicy { UseOriginal, RequirePreservation };

struct CLSEquilibriumAWFReference
{
    double frictionVelocity, traction, matchingVelocity, numerator;
    double production, lambda, minimumLambda, cutoff, logAerodynamicRoughness;
    bool admissible;
    AWFFallbackReason reason;
};

struct SOLVER_UTILS_EXPORT AWFSwitchParameters
{
    // Experimental parameters: no universal terrain calibration is implied.
    double on, off;
    std::array<double, 7> weights;
    AWFVector incomingDirection{1.0, 0.0};
    void Validate() const;
};
struct AWFIndependentSamples
{
    AWFVector matchingVelocity, velocityP, velocityN;
    double heightN, kN, epsilonN, normalVelocityP;
    AWFVector pressureForcing, convectionForcing;
};
struct EquilibriumAWFValues
{
    AWFValues values;
    double requestedWeight, weight, departure;
    std::array<double, 7> sensors;
    AWFFallbackReason reason;
};
SOLVER_UTILS_EXPORT double AWFEquilibriumWeight(double departure,
                                               double on, double off);

/** Dissipation-only correction from the note. EP and EH are values of the same
 * auxiliary distribution as its mean. Momentum, reconstruction, and production
 * always come from Original(). P0 uses current K and ZERO forcing, irrespective
 * of the current flow. An inadmissible cutoff is never clipped.
 */
class SOLVER_UTILS_EXPORT CLSEquilibriumWallFunction
{
public:
    CLSEquilibriumWallFunction(const AWFGeometry &geometry, double kP,
        const RANSNeutralReference &modelReference,
        const AWFParameters &parameters = {});
    const AnalyticalWallFunction &Original() const { return m_original; }
    const CLSEquilibriumAWFReference &Reference() const { return m_reference; }
    EquilibriumAWFValues EvaluateFixedWeight(
        const AWFVector &matchingVelocity, const AWFVector &forcing, double weight,
        AWFFailurePolicy policy = AWFFailurePolicy::UseOriginal,
        bool highReWallRegime = true) const;
    EquilibriumAWFValues Evaluate(
        const AWFIndependentSamples &samples, const AWFSwitchParameters &parameters,
        AWFFailurePolicy policy = AWFFailurePolicy::UseOriginal,
        bool highReWallRegime = true) const;
    double Dissipation(double z, double effectiveWeight) const;

private:
    AnalyticalWallFunction m_original;
    RANSNeutralReference m_model;
    CLSEquilibriumAWFReference m_reference;
};
} // namespace Nektar::SolverUtils
#endif
