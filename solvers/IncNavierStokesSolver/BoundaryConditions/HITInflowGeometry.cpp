///////////////////////////////////////////////////////////////////////////////
// File: HITInflowGeometry.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/BoundaryConditions/HITInflowGeometry.h>
#include <LibUtilities/Foundations/Basis.h>
#include <SpatialDomains/PointGeom.h>
#include <boost/algorithm/string/predicate.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>

namespace Nektar
{
namespace
{
void Require(bool condition, const std::string &message)
{
    if (!condition)
    {
        throw std::invalid_argument(message);
    }
}

// Never let a rank leave setup while its peers enter another collective.
void Agree(const LibUtilities::CommSharedPtr &comm,
           const std::string &localError)
{
    int owner = localError.empty() ? comm->GetSize() : comm->GetRank();
    comm->AllReduce(owner, LibUtilities::ReduceMin);
    if (owner == comm->GetSize())
    {
        return;
    }

    int length = comm->GetRank() == owner ? localError.size() : 0;
    comm->Bcast(length, owner);
    Array<OneD, char> bytes(length);
    if (comm->GetRank() == owner)
    {
        std::copy(localError.begin(), localError.end(), bytes.begin());
    }
    comm->Bcast(bytes, owner);
    NEKERROR(ErrorUtil::efatal,
             "HITInflow geometry: " + std::string(bytes.begin(), bytes.end()));
}

enum RegionKind
{
    Inlet = 1,
    Outlet,
    Periodic
};

std::pair<int, int> CheckConditions(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields, int localRegion,
    int globalRegion, int inletRegion)
{
    using namespace SpatialDomains;
    const auto &first = fields[0]->GetBndConditions()[localRegion];
    if (globalRegion == inletRegion)
    {
        for (int v = 0; v < 3; ++v)
        {
            const auto &condition = fields[v]->GetBndConditions()[localRegion];
            Require(condition->GetBoundaryConditionType() == eDirichlet &&
                        boost::iequals(condition->GetUserDefined(), "HITInflow"),
                    "inlet REF=" + std::to_string(globalRegion) +
                        " requires HITInflow Dirichlet u, v and w");
        }
        const auto &pressure = fields[3]->GetBndConditions()[localRegion];
        Require(pressure->GetBoundaryConditionType() == eNeumann &&
                    boost::iequals(pressure->GetUserDefined(), "H"),
                "HITInflow requires a Neumann pressure condition tagged H");
        return {Inlet, -1};
    }

    for (int v = 0; v < 4; ++v)
    {
        Require(!boost::iequals(
                    fields[v]->GetBndConditions()[localRegion]->GetUserDefined(),
                    "HITInflow"),
                "HITInflow tags occur outside the configured inlet REF");
    }

    if (first->GetBoundaryConditionType() == ePeriodic)
    {
        auto periodic = std::static_pointer_cast<PeriodicBoundaryCondition>(first);
        int partner   = periodic->m_connectedBoundaryRegion;
        for (int v = 0; v < 4; ++v)
        {
            const auto &condition = fields[v]->GetBndConditions()[localRegion];
            Require(condition->GetBoundaryConditionType() == ePeriodic,
                    "periodicity must agree for u, v, w and p");
            Require(!boost::icontains(condition->GetUserDefined(), "Rotated"),
                    "rotational periodicity is unsupported");
            Require(std::static_pointer_cast<PeriodicBoundaryCondition>(condition)
                            ->m_connectedBoundaryRegion ==
                        static_cast<unsigned int>(partner),
                    "periodic partner REF must agree for u, v, w and p");
        }
        return {Periodic, partner};
    }

    for (int v = 0; v < 3; ++v)
    {
        const auto &condition = fields[v]->GetBndConditions()[localRegion];
        Require(condition->GetBoundaryConditionType() == eNeumann ||
                    (condition->GetBoundaryConditionType() == eRobin &&
                     boost::iequals(condition->GetUserDefined(), "HOutflow")),
                "the x outlet requires Neumann or HOutflow Robin velocity");
    }
    const auto &pressure = fields[3]->GetBndConditions()[localRegion];
    Require(pressure->GetBoundaryConditionType() == eDirichlet ||
                (pressure->GetBoundaryConditionType() == eRobin &&
                 boost::iequals(pressure->GetUserDefined(), "HOutflow")),
            "the x outlet requires Dirichlet or HOutflow Robin pressure");
    return {Outlet, -1};
}
} // namespace

std::vector<int> HITInflowGeometry::Validate(
    const LibUtilities::SessionReaderSharedPtr &session,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const HITInflowConfiguration &configuration,
    const SpatialDomains::BoundaryConditions &boundaryConditions)
{
    auto comm = session->GetComm()->GetSpaceComm();
    std::string error;
    std::vector<int> localIds;
    Array<OneD, MultiRegions::ExpListSharedPtr> geometryFields(4);
    bool homogeneous = false;
    bool countGeometry = true;
    int dimensions = 3;
    NekDouble homogeneousLength = 0.0;
    try
    {
        Require(fields.size() == 4,
                "only the u, v, w, p field layout is supported");
        Require(!session->DefinesTag("CreateBndRegions"),
                "CreateBndRegions is unsupported");
        Require(fields[0] &&
                    (fields[0]->GetExpType() == MultiRegions::e3D ||
                     fields[0]->GetExpType() == MultiRegions::e3DH1D),
                "only full 3D and Homogeneous1D fields are supported");
        homogeneous = fields[0]->GetExpType() == MultiRegions::e3DH1D;
        dimensions = homogeneous ? 2 : 3;
        for (int v = 0; v < 4; ++v)
        {
            Require(fields[v] &&
                        fields[v]->GetExpType() == fields[0]->GetExpType(),
                    "velocity and pressure expansion types must agree");
            geometryFields[v] =
                homogeneous ? fields[v]->GetPlane(0) : fields[v];
        }
        if (homogeneous)
        {
            Require(!session->DefinesSolverInfo("HomoStrip"),
                    "HITInflow does not support homogeneous strips");
            homogeneousLength = fields[0]->GetHomoLen();
            Require(std::isfinite(homogeneousLength) && homogeneousLength > 0.0,
                    "the homogeneous period must be finite and positive");
            Require(std::abs(homogeneousLength -
                             configuration.generator.lengths[2]) <=
                        1e-10 * std::max(1.0, homogeneousLength),
                    "synthesis LZ must equal the homogeneous period");
            auto basis = fields[0]->GetHomogeneousBasis();
            const int modes = basis->GetNumModes();
            Require(basis->GetBasisType() == LibUtilities::eFourier &&
                        basis->GetPointsType() ==
                            LibUtilities::eFourierEvenlySpaced &&
                        modes >= 4 && modes % 2 == 0 &&
                        basis->GetNumPoints() == modes,
                    "Homogeneous1D requires an even full Fourier basis with at "
                    "least four modes and matching physical planes");
            Require(configuration.generator.kMax <
                        std::acos(-1.0) * modes / homogeneousLength,
                    "HIT cutoff must lie strictly below the solver's "
                    "homogeneous Nyquist wavenumber");
            const auto planeIds = fields[0]->GetZIDs();
            Require(planeIds.size() > 0,
                    "each homogeneous rank must own physical planes");
            // Every homogeneous rank has the same xy partition. Count its
            // geometric measure once, on the owner of global physical plane 0.
            countGeometry = planeIds[0] == 0;
            for (int v = 0; v < 4; ++v)
            {
                auto otherBasis = fields[v]->GetHomogeneousBasis();
                const auto otherIds = fields[v]->GetZIDs();
                Require(fields[v]->GetHomoLen() == homogeneousLength &&
                            otherBasis->GetBasisType() ==
                                LibUtilities::eFourier &&
                            otherBasis->GetPointsType() ==
                                LibUtilities::eFourierEvenlySpaced &&
                            otherBasis->GetNumModes() == modes &&
                            otherBasis->GetNumPoints() == modes &&
                            otherIds.size() == planeIds.size() &&
                            std::equal(planeIds.begin(), planeIds.end(),
                                       otherIds.begin()),
                        "velocity and pressure homogeneous bases and plane "
                        "ownership must agree");
                Require(geometryFields[v] &&
                            (geometryFields[v]->GetExpSize() == 0 ||
                             geometryFields[v]->GetCoordim(0) == 2),
                        "Homogeneous1D requires an xy mesh with coordinate "
                        "dimension two");
            }
        }

        // DisContField::GenerateBoundaryConditionExpansion visits this map in
        // ascending REF order and includes periodic regions in its arrays.
        for (const auto &region : boundaryConditions.GetBoundaryRegions())
        {
            localIds.push_back(region.first);
        }
        for (int v = 0; v < 4; ++v)
        {
            Require(fields[v]->GetBndConditions().size() == localIds.size() &&
                        fields[v]->GetBndCondExpansions().size() == localIds.size(),
                    "local boundary-array sizes do not match boundary REF IDs");
        }
    }
    catch (const std::exception &e)
    {
        error = e.what();
    }
    Agree(comm, error);

    // A rank can own no inlet or no boundary at all. Build one shared index of
    // all REF IDs before reducing fixed-size statistics on every rank.
    std::set<int> allIdSet(localIds.begin(), localIds.end());
    // CommSerial's gather methods are no-ops, so retain the local list directly
    // for serial setup instead of relying on a serial gather to copy data.
    if (comm->GetSize() > 1)
    {
        Array<OneD, int> localCount(1, localIds.size()), counts(comm->GetSize());
        comm->AllGather(localCount, counts);
        Array<OneD, int> offsets(comm->GetSize(), 0);
        for (int r = 1; r < comm->GetSize(); ++r)
        {
            offsets[r] = offsets[r - 1] + counts[r - 1];
        }
        int total = offsets[comm->GetSize() - 1] + counts[comm->GetSize() - 1];
        Array<OneD, int> localArray(localIds.size()), allArray(total);
        std::copy(localIds.begin(), localIds.end(), localArray.begin());
        comm->AllGatherv(localArray, allArray, counts, offsets);
        allIdSet.insert(allArray.begin(), allArray.end());
    }
    std::vector<int> globalIds(allIdSet.begin(), allIdSet.end());
    std::map<int, int> index;
    for (int i = 0; i < globalIds.size(); ++i)
    {
        index[globalIds[i]] = i;
    }

    const int n = globalIds.size();
    const NekDouble large = std::numeric_limits<NekDouble>::max();
    Array<OneD, NekDouble> domainMin(3, large), domainMax(3, -large);
    // First three entries are coordinate extrema, next three normal extrema.
    Array<OneD, NekDouble> regionMin(6 * n, large), regionMax(6 * n, -large);
    Array<OneD, NekDouble> area(n, 0.0);
    Array<OneD, int> points(n, 0);
    Array<OneD, int> metadataMin(2 * n, std::numeric_limits<int>::max());
    Array<OneD, int> metadataMax(2 * n, std::numeric_limits<int>::min());
    try
    {
        Require(n == 2 * dimensions,
                homogeneous
                    ? "Homogeneous1D requires four rectangular xy boundary "
                      "regions and implicit periodicity in z"
                    : "full 3D requires one rectangular boundary region on "
                      "each of the six box sides");
        Require(index.count(configuration.regionId) == 1,
                "configured inlet REF is absent from the mesh");

        Array<OneD, Array<OneD, NekDouble>> coordinates(3);
        for (int d = 0; d < 3; ++d)
        {
            coordinates[d] =
                Array<OneD, NekDouble>(geometryFields[0]->GetTotPoints(), 0.0);
        }
        geometryFields[0]->GetCoords(coordinates[0], coordinates[1],
                                    coordinates[2]);
        for (int d = 0; d < dimensions; ++d)
        {
            for (NekDouble value : coordinates[d])
            {
                Require(std::isfinite(value), "non-finite mesh coordinate");
                domainMin[d] = std::min(domainMin[d], value);
                domainMax[d] = std::max(domainMax[d], value);
            }
        }
        // Collapsed-coordinate simplex quadrature need not contain every
        // vertex. Include vertices so a valid tetrahedral box has exact spans.
        for (int e = 0; e < geometryFields[0]->GetExpSize(); ++e)
        {
            auto geometry = geometryFields[0]->GetExp(e)->GetGeom();
            for (int v = 0; v < geometry->GetNumVerts(); ++v)
            {
                std::array<NekDouble, 3> point;
                geometry->GetVertex(v)->GetCoords(point[0], point[1], point[2]);
                for (int d = 0; d < dimensions; ++d)
                {
                    Require(std::isfinite(point[d]), "non-finite mesh vertex");
                    domainMin[d] = std::min(domainMin[d], point[d]);
                    domainMax[d] = std::max(domainMax[d], point[d]);
                }
            }
        }
        if (homogeneous)
        {
            // Fourier quadrature excludes the duplicated z=LZ endpoint. Its
            // geometric extent is the implicit period, not sampled z extrema.
            domainMin[2] = 0.0;
            domainMax[2] = homogeneousLength;
        }

        for (int local = 0; local < localIds.size(); ++local)
        {
            int global = index.at(localIds[local]);
            auto kind = CheckConditions(fields, local, localIds[local],
                                        configuration.regionId);
            metadataMin[2 * global] = metadataMax[2 * global] = kind.first;
            metadataMin[2 * global + 1] = metadataMax[2 * global + 1] = kind.second;

            auto boundary = geometryFields[0]->GetBndCondExpansions()[local];
            Require(boundary && boundary->GetTotPoints() > 0,
                    "empty local boundary expansion");
            for (int v = 1; v < 4; ++v)
            {
                auto other = geometryFields[v]->GetBndCondExpansions()[local];
                Require(other && other->GetExpSize() == boundary->GetExpSize(),
                        "boundary geometry differs between velocity and pressure");
                for (int e = 0; e < boundary->GetExpSize(); ++e)
                {
                    Require(other->GetExp(e)->GetGeom()->GetGlobalID() ==
                                boundary->GetExp(e)->GetGeom()->GetGlobalID(),
                            "boundary geometry order differs between fields");
                }
            }

            for (int d = 0; d < 3; ++d)
            {
                coordinates[d] =
                    Array<OneD, NekDouble>(boundary->GetTotPoints(), 0.0);
            }
            boundary->GetCoords(coordinates[0], coordinates[1], coordinates[2]);
            Array<OneD, Array<OneD, NekDouble>> normals;
            geometryFields[0]->GetBoundaryNormals(local, normals);
            for (int d = 0; d < dimensions; ++d)
            {
                for (int q = 0; q < boundary->GetTotPoints(); ++q)
                {
                    Require(std::isfinite(coordinates[d][q]) &&
                                std::isfinite(normals[d][q]),
                            "non-finite boundary coordinate or normal");
                    regionMin[6 * global + d] = std::min(
                        regionMin[6 * global + d], coordinates[d][q]);
                    regionMax[6 * global + d] = std::max(
                        regionMax[6 * global + d], coordinates[d][q]);
                    regionMin[6 * global + d + 3] = std::min(
                        regionMin[6 * global + d + 3], normals[d][q]);
                    regionMax[6 * global + d + 3] = std::max(
                        regionMax[6 * global + d + 3], normals[d][q]);
                }
            }
            if (homogeneous)
            {
                regionMin[6 * global + 2] = 0.0;
                regionMax[6 * global + 2] = homogeneousLength;
                regionMin[6 * global + 5] = 0.0;
                regionMax[6 * global + 5] = 0.0;
            }
            points[global] = boundary->GetTotPoints();

            for (int e = 0; e < boundary->GetExpSize(); ++e)
            {
                auto geometry = boundary->GetExp(e)->GetGeom();
                for (int v = 0; v < geometry->GetNumVerts(); ++v)
                {
                    std::array<NekDouble, 3> point;
                    geometry->GetVertex(v)->GetCoords(point[0], point[1], point[2]);
                    for (int d = 0; d < dimensions; ++d)
                    {
                        Require(std::isfinite(point[d]),
                                "non-finite boundary vertex");
                        regionMin[6 * global + d] =
                            std::min(regionMin[6 * global + d], point[d]);
                        regionMax[6 * global + d] =
                            std::max(regionMax[6 * global + d], point[d]);
                    }
                }
            }

            // ExpList::Integral is collective on the boundary communicator.
            // Integrate individual faces locally, then reduce exactly once.
            Array<OneD, NekDouble> ones(boundary->GetTotPoints(), 1.0);
            for (int e = 0; countGeometry && e < boundary->GetExpSize(); ++e)
            {
                area[global] +=
                    (homogeneous ? homogeneousLength : 1.0) *
                    boundary->GetExp(e)->Integral(
                        ones + boundary->GetPhys_Offset(e));
            }
        }
    }
    catch (const std::exception &e)
    {
        error = e.what();
    }
    Agree(comm, error);
    comm->AllReduce(domainMin, LibUtilities::ReduceMin);
    comm->AllReduce(domainMax, LibUtilities::ReduceMax);
    comm->AllReduce(regionMin, LibUtilities::ReduceMin);
    comm->AllReduce(regionMax, LibUtilities::ReduceMax);
    comm->AllReduce(area, LibUtilities::ReduceSum);
    comm->AllReduce(points, LibUtilities::ReduceSum);
    comm->AllReduce(metadataMin, LibUtilities::ReduceMin);
    comm->AllReduce(metadataMax, LibUtilities::ReduceMax);

    try
    {
        std::array<NekDouble, 3> lengths;
        NekDouble scale = 1.0;
        for (int d = 0; d < 3; ++d)
        {
            lengths[d] = domainMax[d] - domainMin[d];
            Require(std::isfinite(lengths[d]) && lengths[d] > 0.0,
                    "the mesh must have a finite positive extent in each axis");
            scale = std::max(scale, std::abs(domainMin[d]));
            scale = std::max(scale, std::abs(domainMax[d]));
        }
        const NekDouble tolerance =
            1e-9 * *std::max_element(lengths.begin(), lengths.end()) +
            128.0 * std::numeric_limits<NekDouble>::epsilon() * scale;
        const auto near = [tolerance](NekDouble a, NekDouble b)
        { return std::abs(a - b) <= tolerance; };
        Require(near(configuration.generator.origin[0], domainMin[0]),
                "ORIGIN.X must equal the upstream x boundary");
        for (int d = 1; d < 3; ++d)
        {
            Require(near(configuration.generator.lengths[d], lengths[d]),
                    "synthesis LY and LZ must equal the mesh transverse periods");
        }

        std::array<int, 6> sides;
        sides.fill(-1);
        for (int i = 0; i < n; ++i)
        {
            const std::string label = "REF=" + std::to_string(globalIds[i]);
            Require(points[i] > 0, label + " has no boundary points");
            Require(metadataMin[2 * i] == metadataMax[2 * i] &&
                        metadataMin[2 * i + 1] == metadataMax[2 * i + 1],
                    label + " boundary types or partners disagree across ranks");
            int side = -1;
            for (int d = 0; d < dimensions; ++d)
            {
                if (!near(regionMin[6 * i + d], regionMax[6 * i + d]))
                {
                    continue;
                }
                for (int upper = 0; upper < 2; ++upper)
                {
                    NekDouble coordinate = upper ? domainMax[d] : domainMin[d];
                    if (near(regionMin[6 * i + d], coordinate) &&
                        near(regionMax[6 * i + d], coordinate))
                    {
                        Require(side < 0, label + " has degenerate geometry");
                        side = 2 * d + upper;
                    }
                }
            }
            Require(side >= 0, label + " must be a planar axis-aligned box side");
            Require(sides[side] < 0, "each box side must have one boundary REF");
            sides[side] = i;
            int axis = side / 2;
            for (int d = 0; d < 3; ++d)
            {
                NekDouble expectedNormal = d == axis ? (side % 2 ? 1.0 : -1.0) : 0.0;
                Require(std::abs(regionMin[6 * i + d + 3] - expectedNormal) < 1e-8 &&
                            std::abs(regionMax[6 * i + d + 3] - expectedNormal) < 1e-8,
                        label + " has an unsupported exterior normal");
                if (d != axis)
                {
                    Require(near(regionMin[6 * i + d], domainMin[d]) &&
                                near(regionMax[6 * i + d], domainMax[d]),
                            label + " must span the complete box side");
                }
            }
            NekDouble expectedArea = lengths[(axis + 1) % 3] *
                                     lengths[(axis + 2) % 3];
            Require(std::isfinite(area[i]) &&
                        std::abs(area[i] - expectedArea) <= 1e-7 * expectedArea,
                    label + " boundary area does not cover the box side");
            int expectedKind = side == 0 ? Inlet : side == 1 ? Outlet : Periodic;
            Require(metadataMin[2 * i] == expectedKind,
                    label + " must be the x- inlet, x+ outlet, or a transverse "
                            "periodic boundary according to its location");
            Require(side != 0 || globalIds[i] == configuration.regionId,
                    "configured HITInflow REF must be the x- box side");
        }

        for (int side = 0; side < 2 * dimensions; ++side)
        {
            Require(sides[side] >= 0, "a rectangular box side is missing");
            if (side >= 2)
            {
                int i = sides[side];
                int opposite = sides[side ^ 1];
                Require(metadataMin[2 * i + 1] == globalIds[opposite] &&
                            metadataMin[2 * opposite + 1] == globalIds[i],
                        "transverse periodic partners must be reciprocal "
                        "opposite sides on the same axis");
            }
        }
    }
    catch (const std::exception &e)
    {
        error = e.what();
    }
    Agree(comm, error);
    return localIds;
}
} // namespace Nektar
