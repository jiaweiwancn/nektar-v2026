///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_RANS_WALL_GLOBAL_FIELD_H
#define NEKTAR_RANS_WALL_GLOBAL_FIELD_H

#include <IncNavierStokesSolver/EquationSystems/RANSWallFieldLayout.h>
#include <array>
#include <functional>
#include <map>
#include <memory>
#include <tuple>
#include <vector>

namespace Nektar
{
using RANSWallCoordinate = std::array<double,3>;
using RANSWallFunction = std::function<double(const RANSWallCoordinate &)>;

/** A row of a serial native CG space in physical Cartesian coordinates.
 * The constant vector is the native coefficient representation of one on
 * this element, including local/global signs. Value and Gradient are centred
 * on this constant, so constant fields evaluate and differentiate exactly.
 */
struct RANSWallGlobalPoint
{
    using Vector = RANSWallAlgebra::Vector;
    RANSWallCoordinate coordinate{};
    int geometry = -1, element = -1;
    std::vector<int> ids;
    Vector phi;
    std::array<Vector,3> derivative;
    Vector constant;
    double weight = 0;
    bool columnAverage = false; // normalized physical functional, not a point
    double Value(const Vector &state, int offset=0) const;
    std::array<double,3> Gradient(const Vector &state, int offset=0) const;
};

/** Native tensor CG field on canonical affine quads or hexahedra.
 * The caller supplies a serial session on each MPI rank. This class owns no
 * MPI collectives; the patch assigns physical element/quadrature ownership.
 * Physical axes are x,z for dimension two and x,y,z for dimension three.
 */
class RANSWallGlobalField
{
public:
    using Coordinate = RANSWallCoordinate;
    using Vector = RANSWallAlgebra::Vector;
    using WallMatrix = RANSWallAlgebra::DenseMatrix;
    struct Element
    {
        Coordinate lower{}, upper{};
        int geometry = -1, degree = 0;
    };
    using BoundaryPredicate = std::function<bool(int,int,const Element &)>;

    RANSWallGlobalField(const LibUtilities::SessionReaderSharedPtr &serialSession,
                        const MultiRegions::ExpListSharedPtr &full,
                        const std::vector<unsigned int> &elements, int dimension);
    const MultiRegions::ContFieldSharedPtr &Field() const { return m_field; }
    const std::vector<Element> &Elements() const { return m_elements; }
    int Dimension() const { return m_dimension; }
    int Size() const;
    RANSWallGlobalPoint Point(const Coordinate &coordinate, int side=0) const;
    int LocateElement(const Coordinate &coordinate,int side=0) const;
    RANSWallGlobalPoint ElementPoint(int element, const Coordinate &coordinate) const;
    /** d/dz of each Cartesian basis-gradient row on the same native element.
     * Used by the opt-in R3 action when a normalization-ray point moves.
     * Contracting with a field uses the same constant-centering as Gradient.
     */
    std::array<Vector,3> RadialGradientDerivative(const RANSWallGlobalPoint &point) const;
    std::vector<RANSWallGlobalPoint> Quadrature(int order=0) const;
    WallMatrix Mass() const;
    Vector Project(const RANSWallFunction &function) const;
    std::map<int,double> Boundary(const RANSWallFunction &function,
                                  const BoundaryPredicate &predicate) const;
    void Store(const Vector &state, int offset=0) const;

    /** Constant-horizontal physical L2 projection with value/derivative KKT
     * constraints {height,side,value}. side=0 selects a value; -1/+1 select
     * a physical z derivative approached from below/above. The returned full
     * native vector has only tangential endpoint modes; all higher tangential
     * modes are exactly zero. Radial interval partitions must conform across
     * columns. Variable horizontal areas enter the physical 1D measure.
     */
    Vector RadialProject(const std::function<double(double)> &function,
        const std::vector<std::array<double,3>> &constraints,
        const std::vector<double> &physicalBreakpoints={}) const;

    /** Independent diagnostic through full native BwdTrans/PhysDeriv and
     * standard physical interpolation, not the fast tensor evaluation rows.
     * Returns value, dx,dy,dz in physical Cartesian axes.
     */
    /** Explicit element ownership for independent derivative replay on a
     * vertical face or H/corner intersection. No coordinate side is overloaded.
     */
    std::array<double,4> NativeEvaluateElement(const Vector &state,int element,
        const Coordinate &coordinate,int offset=0) const;
    std::array<double,4> NativeEvaluate(const Vector &state,
                                       const Coordinate &coordinate,
                                       int offset=0, int side=0) const;

private:
    struct BasisData;
    MultiRegions::ContFieldSharedPtr m_field;
    int m_dimension;
    std::array<int,3> m_axes{};
    std::vector<Element> m_elements;
    std::vector<std::array<std::shared_ptr<BasisData>,3>> m_elementBasis;
    std::map<std::tuple<int,int,int>,std::shared_ptr<BasisData>> m_basis;
};
}
#endif
