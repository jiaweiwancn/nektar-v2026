///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SOLVERUTILS_RANS_WALL_NORMAL_MAP_H
#define NEKTAR_SOLVERUTILS_RANS_WALL_NORMAL_MAP_H

#include <SolverUtils/Turbulence/RANSModels.h>
#include <string>

namespace Nektar::SolverUtils
{
/** Analytic geometry over a canonical native tensor chart (q1,q2,rho).
 * Physical vectors always use Cartesian (x,y,z). In 2D, q2=y=0 and the
 * geometry is extruded in y; the 3x3 metrics retain that unit span direction.
 * H(q) is physical matching distance and height is its reference/chart value.
 * No CLS/AWF constitutive law is changed by this coordinate mapping.
 */
struct SOLVER_UTILS_EXPORT RANSWallNormalMap
{
    enum class Kind { Flat, Plane, Cylinder, Gaussian, CosineHill };
    int dimension = 2;
    Kind kind = Kind::Flat;
    double height = 2, heightSlopeX = 0, heightSlopeY = 0;
    double heightAmplitude = 0, heightWaveX = 1, heightWaveY = 1;
    double slopeX = 0, slopeY = 0, radius = 20;
    double amplitude = .02, lengthX = 1, lengthY = 1;

    struct Sample
    {
        RANSVector chart{}, point{}, normal{}, tangent1{}, tangent2{};
        // jacobian[a][b]=dx_a/dchart_b; inverse[b][a]=dchart_b/dx_a.
        // B=grad_Gamma(n), Br=B*(I+rB)^-1; both act tangentially.
        RANSTensor jacobian{}, inverse{}, B{}, Br{};
        double determinant = 0, baseArea = 0, distance = 0, height = 0;
        double scale = 0, J = 0, kappa1 = 0, kappa2 = 0;
    };

    /** At an exact cosine foot, side=-1/+1 selects the q1 left/right limit;
     * side=0 selects the flat-side curvature. Nonzero side elsewhere is an
     * error. Explicit foot patches must select their own one-sided metric.
     */
    Sample At(const RANSVector &chart, int footSide = 0) const;

    /** Unique normal projection is required. Cylinder q1 uses the principal
     * atan2 branch [-pi*radius,pi*radius]; unwrap into a certified nonwrapping
     * chart interval when it crosses this branch cut. No domain is cached.
     */
    RANSVector Inverse(const RANSVector &point) const;

    // Oriented physical face area vector per unit of the other chart axes.
    RANSVector Cofactor(const RANSVector &chart, int axis, int side = 1,
                        int footSide = 0) const;
    Sample Matching(const RANSVector &chart) const;
    RANSVector MatchingArea(const RANSVector &chart, int footSide = 0) const;
    bool Identity() const;
    bool PlanarConstantHeight() const;

    /** Sufficient whole-tube certificate, including the outer rho range.
     * Graphs use a global C1,1 quadratic-remainder bound; the outward cylinder
     * uses its nonwrapping arc. This assumes an isolated wall: clearance from
     * other walls is a separate requirement. C1 cosine feet remain explicit
     * patch interfaces. Physical roughness/rK/rm checks belong to the caller.
     */
    void Validate(const RANSVector &lower, const RANSVector &upper,
                  double minimumStretch = .2) const;
    std::string Signature() const;
};
} // namespace Nektar::SolverUtils
#endif
