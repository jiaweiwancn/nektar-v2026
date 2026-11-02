///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_RANS_WALL_INTERFACE_H
#define NEKTAR_RANS_WALL_INTERFACE_H
#include <IncNavierStokesSolver/EquationSystems/RANSWallPatch.h>
namespace Nektar
{
/** Matching flat native traces. This is the conforming case: stage velocity
 * and k use the full-domain native CG injection, while pressure stays broken.
 * A common numerical flux has opposite physical dual loads in the two spaces.
 */
class RANSWallInterface
{
public:
    using Vector=RANSWallAlgebra::Vector;
    using WallMatrix=RANSWallAlgebra::DenseMatrix;
    explicit RANSWallInterface(const RANSWallFieldLayout &layout,int order=12);
    const WallMatrix &TraceConstraint() const { return m_traceConstraint; }
    const std::vector<int> &TraceIds() const { return m_traceIds; }
    const std::vector<RANSWallPoint> &Points() const { return m_points; }
    // Pair of native subdomain dual loads from ONE supplied physical flux.
    std::array<Vector,2> Loads(const RANSWallAlgebra::Function &flux) const;
    Vector Merge(const std::array<Vector,2> &loads) const;
    struct PressureTest
    {
        WallMatrix divergence,system;
        Vector rhs,solution;
        std::vector<int> freeVelocity;
        double velocityMaximum,globalMean,pressureJump,tractionMismatch;
    };
    PressureTest ReferencePressure(const SolverUtils::RANSNeutralReference &model,double k) const;
private:
    const RANSWallFieldLayout &m_layout;
    std::vector<RANSWallPoint> m_points;
    std::vector<int> m_traceIds;
    WallMatrix m_traceConstraint;
    std::array<std::vector<std::pair<int,double>>,2> m_injection;
};
}
#endif
