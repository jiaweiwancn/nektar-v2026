///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SOLVERUTILS_RANS_WALL_GEOMETRY_H
#define NEKTAR_SOLVERUTILS_RANS_WALL_GEOMETRY_H

#include <SolverUtils/Turbulence/AnalyticalWallFunction.h>

namespace Nektar::SolverUtils
{
enum class RANSTerrain { Flat, Plane, CosineHill };

/** Physical x/span/z frame. A 2D solver maps its (x,y) to physical (x,z).
 * The terrain is a graph, extruded in span; eta is physical normal distance.
 * Curvature is signed: dt/ds = curvature*n, n points into the fluid.
 */
struct SOLVER_UTILS_EXPORT RANSWallFrame
{
    double x, z, distance, curvature;
    RANSVector tangent, normal;
    double Jacobian(double eta) const;
    double Volume(double height) const;
};

struct SOLVER_UTILS_EXPORT RANSWallGeometry
{
    RANSTerrain terrain = RANSTerrain::Flat;
    double origin = 0, centre = 0, slope = 0, hillHeight = .04, hillHalfLength = .1;
    void Validate() const;
    RANSWallFrame At(double footX) const;
    RANSWallFrame Closest(double x, double z) const;
    RANSVector Point(double footX, double eta, double span = 0) const;
    void CheckStrip(const RANSWallFrame &frame, double height, double pointN,
                    double maximumCurvatureRatio = .25) const;
};

struct SOLVER_UTILS_EXPORT AWFSpatialTransition
{
    bool enabled = false;
    double begin = 0, end = 1;
    void Validate() const;
    double Weight(double footX) const;
};

struct AWFMetricValues
{
    double meanProduction, meanDissipation;
    AWFVector integratedShear;
};

/** Integrate the original auxiliary fields with J=1-curvature*eta. The
 * correction, when nonzero, is the same capped distribution used at P and H.
 * This changes geometric integration, not the local AWF momentum model.
 */
SOLVER_UTILS_EXPORT AWFMetricValues IntegrateAWFMetric(
    const AnalyticalWallFunction &wall, const AWFVector &velocity,
    const AWFVector &forcing, double curvature, double weight = 0,
    double correctedNumerator = 0, double correctedCutoff = 0);
}
#endif
