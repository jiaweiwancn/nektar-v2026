///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_RANS_WALL_FIELD_LAYOUT_H
#define NEKTAR_RANS_WALL_FIELD_LAYOUT_H
#include <MultiRegions/ContField.h>
#include <array>
#include <functional>
#include <map>
#include <vector>

namespace Nektar
{
namespace RANSWallAlgebra
{
using Vector = std::vector<double>;
using Matrix = std::vector<Vector>;
using Function = std::function<double(double, double)>;
Vector Solve(const Matrix &matrix, const Vector &rhs);
Vector Multiply(const Matrix &matrix, const Vector &value);
double Maximum(const Vector &value);
void Require(bool valid, const std::string &message);
}

/** One physical evaluation row in a native continuous expansion. Derivatives
 * are physical mesh (x,y)=(x,z) derivatives, evaluated by LocalRegions. */
struct RANSWallPoint
{
    double x = 0, r = 0, weight = 0;
    int element = -1, geometry = -1;
    std::vector<int> ids;
    RANSWallAlgebra::Vector phi, dx, dr;
    double Value(const RANSWallAlgebra::Vector &field, int offset = 0) const;
    std::array<double, 2> Gradient(const RANSWallAlgebra::Vector &field, int offset = 0) const;
};

/** Native ExpList/ContField/AssemblyMapCG with independent physical quadrature.
 * The M3 backend accepts serial affine axis-aligned 2D quads. Unsupported
 * geometry/partitioning is rejected; it is not silently treated as flat.
 * It never asks a full-domain epsilon field to mask out unwanted PDE rows. */
class RANSWallNativeField
{
public:
    struct Element
    {
        double left, right, bottom, top;
        int geometry, degree;
        std::vector<RANSWallAlgebra::Vector> derivativeX, derivativeR;
    };
    RANSWallNativeField(const LibUtilities::SessionReaderSharedPtr &session,
        const MultiRegions::ExpListSharedPtr &full, const std::vector<unsigned int> &elements);
    int Size() const;
    const MultiRegions::ContFieldSharedPtr &Field() const { return m_field; }
    const std::vector<Element> &Elements() const { return m_elements; }
    RANSWallPoint Point(double x, double r, int side = 0) const;
    RANSWallPoint ElementPoint(int element, double x, double r) const;
    std::vector<RANSWallPoint> Quadrature(int order = 0) const;
    RANSWallAlgebra::Matrix Mass() const;
    RANSWallAlgebra::Vector Project(const RANSWallAlgebra::Function &function) const;
    std::map<int, double> Boundary(const RANSWallAlgebra::Function &function,
        const std::function<bool(int, const Element &)> &select) const;
    void Store(const RANSWallAlgebra::Vector &state, int offset = 0) const;
    RANSWallAlgebra::Vector NativeLoad(const RANSWallAlgebra::Function &function) const;
private:
    MultiRegions::ContFieldSharedPtr m_field;
    std::vector<Element> m_elements;
};

/** Separate native CG subdomain spaces with common velocity/k trace injection.
 * U,W are full stage blocks; their inner free modes can be condensed. Pressure
 * has two CG blocks and only one global gauge. Epsilon has outer rows only. */
class RANSWallFieldLayout
{
public:
    enum Block { U, W, K, Epsilon, InnerPressure, OuterPressure, Gauge, Blocks };
    RANSWallFieldLayout(const LibUtilities::SessionReaderSharedPtr &session,
        const SpatialDomains::MeshGraphSharedPtr &graph, double height);
    std::shared_ptr<RANSWallNativeField> velocity, innerVelocity, outerVelocity;
    std::shared_ptr<RANSWallNativeField> k, epsilon, innerPressure, outerPressure;
    std::array<int, Blocks> offset{}, size{};
    int Size() const { return offset[Gauge]+1; }
    double Height() const { return m_height; }
    double Left() const { return m_left; }
    double Right() const { return m_right; }
    double Top() const { return m_top; }
    const RANSWallNativeField &Field(Block block) const;
private:
    double m_height, m_left, m_right, m_top;
};
}
#endif
