///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_RANS_PHYSICAL_POINT_SAMPLER_H
#define NEKTAR_RANS_PHYSICAL_POINT_SAMPLER_H
#include <MultiRegions/ExpList.h>
#include <LibUtilities/Communication/Comm.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>

namespace Nektar
{
/** Read-only native interpolation in volume or embedded boundary expansions.
 * Every rank supplies identical query coordinates in mesh x/y/z axes. Global
 * geometry ID, then rank, selects a deterministic one-sided trace at shared
 * points. No nearest quadrature-point or nearest-element extrapolation.
 */
class RANSPhysicalPointSampler
{
public:
    using Coordinate = std::array<double,3>;
    using Fields = Array<OneD,Array<OneD,NekDouble>>;
    RANSPhysicalPointSampler(const std::vector<LocalRegions::ExpansionSharedPtr> &expansions,
        const std::vector<int> &offsets, const LibUtilities::CommSharedPtr &comm,
        const std::vector<Coordinate> &points, int coordinateDimension)
        : m_expansions(expansions),m_offsets(offsets),m_comm(comm),m_points(points),
          m_elements(points.size(),-1),m_local(points.size()),
          m_ids(points.size(),std::numeric_limits<int>::max()),
          m_owners(points.size(),std::numeric_limits<int>::max()),m_residual(points.size(),0.)
    {
        ASSERTL0(expansions.size() == offsets.size(),"Point sampler offset mismatch");
        for (const auto &exp : expansions)
        {
            const auto shape=exp->GetGeom()->GetShapeType();
            ASSERTL0(shape==LibUtilities::eSegment || shape==LibUtilities::eQuadrilateral ||
                     shape==LibUtilities::eHexahedron,
                     "Physical wall observation supports segment/quad/hex tensor-product geometry");
        }
        for (int q = 0; q < points.size(); ++q)
        {
            ASSERTL0(std::all_of(points[q].begin(),points[q].end(),[](double v){return std::isfinite(v);}),
                     "Nonfinite physical observation point");
            for (int e = 0; e < expansions.size(); ++e)
            {
                auto geom = expansions[e]->GetGeom();
                const auto box = geom->GetBoundingBox();
                bool candidate = true;
                for (int i = 0; i < coordinateDimension; ++i)
                    candidate &= points[q][i] >= box[i]-1e-8 && points[q][i] <= box[i+3]+1e-8;
                if (!candidate || geom->GetGlobalID() >= m_ids[q]) continue;
                Array<OneD,NekDouble> local(coordinateDimension,0.);
                double residual;
                if (Inverse(geom,points[q],coordinateDimension,local,residual))
                {
                    m_ids[q] = geom->GetGlobalID();m_elements[q] = e;
                    m_local[q] = local;m_residual[q] = residual;
                }
            }
        }
        m_comm->AllReduce(m_ids,LibUtilities::ReduceMin);
        for (int q = 0; q < points.size(); ++q)
            if (m_elements[q] >= 0 && expansions[m_elements[q]]->GetGeom()->GetGlobalID() == m_ids[q])
                m_owners[q] = m_comm->GetRank();
        m_comm->AllReduce(m_owners,LibUtilities::ReduceMin);
        for (int q = 0; q < points.size(); ++q)
            if (m_owners[q] != m_comm->GetRank()) { m_elements[q] = -1;m_residual[q] = 0.; }
        m_comm->AllReduce(m_residual,LibUtilities::ReduceMax);
    }

    static RANSPhysicalPointSampler Volume(const MultiRegions::ExpListSharedPtr &field,
        const std::vector<Coordinate> &points,int dimension)
    {
        std::vector<LocalRegions::ExpansionSharedPtr> exps;
        std::vector<int> offsets;
        for (int e = 0; e < field->GetExpSize(); ++e)
        { exps.push_back(field->GetExp(e));offsets.push_back(field->GetPhys_Offset(e)); }
        return {exps,offsets,field->GetComm(),points,dimension};
    }
    Fields Sample(const Fields &values) const
    {
        Fields result(values.size());
        for (int f = 0; f < values.size(); ++f)
        {
            result[f] = Array<OneD,NekDouble>(m_points.size(),0.);
            for (int q = 0; q < m_points.size(); ++q)
                if (m_elements[q] >= 0)
                    result[f][q] = m_expansions[m_elements[q]]->StdPhysEvaluate(
                        m_local[q],values[f]+m_offsets[m_elements[q]]);
            m_comm->AllReduce(result[f],LibUtilities::ReduceSum);
        }
        return result;
    }
    bool Found(int q) const { return m_ids[q] != std::numeric_limits<int>::max(); }
    int Geometry(int q) const { return Found(q) ? m_ids[q] : -1; }
    int Owner(int q) const { return Found(q) ? m_owners[q] : -1; }
    double Residual(int q) const { return m_residual[q]; }

private:
    static bool Inverse(SpatialDomains::Geometry *geom,const Coordinate &point,int dim,
        Array<OneD,NekDouble> &local,double &residual)
    {
        Array<OneD,NekDouble> target(dim);
        for (int i = 0; i < dim; ++i) target[i] = point[i];
        const int natural = geom->GetShapeDim();
        // Nektar's Geometry1D::GetLocCoords only implements regular segments.
        // Curved embedded edges use the centre as a Newton starting point;
        // the residual below always tests all ambient coordinates.
        if (natural > 1) geom->GetLocCoords(target,local);
        if (std::any_of(local.begin(),local.end(),[](double v){return !std::isfinite(v);}))
            Vmath::Zero(local.size(),local,1);
        for (int iteration = 0; iteration < 25; ++iteration)
        {
            double r[3]{},j[3][3]{},a[3][4]{};residual = 0;
            for (int i = 0; i < dim; ++i)
            {
                r[i] = point[i]-geom->GetCoord(i,local);
                if (!std::isfinite(r[i])) return false;
                residual = std::max(residual,std::abs(r[i]));
            }
            if (residual < 5e-12)
            {
                for (int k = 0; k < natural; ++k) if (std::abs(local[k]) > 1+1e-8) return false;
                return true;
            }
            for (int k = 0; k < natural; ++k)
            {
                const double saved = local[k],h = 1e-5;
                local[k] = saved+h;
                for (int i = 0; i < dim; ++i) j[i][k] = geom->GetCoord(i,local);
                local[k] = saved-h;
                for (int i = 0; i < dim; ++i) j[i][k] = (j[i][k]-geom->GetCoord(i,local))/(2*h);
                local[k] = saved;
            }
            for (int k = 0; k < natural; ++k)
            {
                for (int i = 0; i < dim; ++i) a[k][natural] += j[i][k]*r[i];
                for (int l = 0; l < natural; ++l)
                    for (int i = 0; i < dim; ++i) a[k][l] += j[i][k]*j[i][l];
            }
            for (int k = 0; k < natural; ++k)
            {
                int pivot = k;
                for (int l = k+1; l < natural; ++l) if (std::abs(a[l][k]) > std::abs(a[pivot][k])) pivot = l;
                if (std::abs(a[pivot][k]) < 1e-24) return false;
                for (int l = k; l <= natural; ++l) std::swap(a[k][l],a[pivot][l]);
                const double diagonal = a[k][k];
                for (int l = k; l <= natural; ++l) a[k][l] /= diagonal;
                for (int row = 0; row < natural; ++row) if (row != k)
                {
                    const double factor = a[row][k];
                    for (int l = k; l <= natural; ++l) a[row][l] -= factor*a[k][l];
                }
            }
            for (int k = 0; k < natural; ++k)
            {
                local[k] += a[k][natural];
                if (!std::isfinite(local[k]) || std::abs(local[k]) > 3) return false;
            }
        }
        return false;
    }
    std::vector<LocalRegions::ExpansionSharedPtr> m_expansions;
    std::vector<int> m_offsets;
    LibUtilities::CommSharedPtr m_comm;
    std::vector<Coordinate> m_points;
    std::vector<int> m_elements;
    std::vector<Array<OneD,NekDouble>> m_local;
    Array<OneD,int> m_ids,m_owners;
    Array<OneD,NekDouble> m_residual;
};
}
#endif
