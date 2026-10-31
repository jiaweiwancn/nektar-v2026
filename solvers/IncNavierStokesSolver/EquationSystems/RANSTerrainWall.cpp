///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSTerrainWall.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/ContField.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace Nektar
{
using namespace SolverUtils;
using MultiRegions::DirCartesianMap;
namespace
{
double Dot(const RANSVector &a, const RANSVector &b)
{
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
bool RefineCoordinates(SpatialDomains::Geometry *geometry,
    const Array<OneD,NekDouble> &point, Array<OneD,NekDouble> &local,
    double &residual, double insideTolerance)
{
    geometry->GetLocCoords(point,local);
    const int dim = point.size();
    for (int iteration = 0; iteration < 20; ++iteration)
    {
        double a[3][4]{}; residual = 0;
        for (int i = 0; i < dim; ++i)
        {
            a[i][dim] = point[i]-geometry->GetCoord(i,local);
            residual = std::max(residual,std::abs(a[i][dim]));
        }
        if (residual < 2e-13)
        {
            for (double v : local) if (std::abs(v) > 1+insideTolerance) return false;
            return true;
        }
        for (int j = 0; j < dim; ++j)
        {
            const double saved = local[j], step = 1e-5;
            local[j] = saved+step;
            for (int i = 0; i < dim; ++i) a[i][j] = geometry->GetCoord(i,local);
            local[j] = saved-step;
            for (int i = 0; i < dim; ++i) a[i][j] = (a[i][j]-geometry->GetCoord(i,local))/(2*step);
            local[j] = saved;
        }
        for (int j = 0; j < dim; ++j)
        {
            int pivot = j;
            for (int i = j+1; i < dim; ++i) if (std::abs(a[i][j]) > std::abs(a[pivot][j])) pivot = i;
            if (std::abs(a[pivot][j]) < 1e-18) return false;
            for (int k = j; k <= dim; ++k) std::swap(a[pivot][k],a[j][k]);
            const double d = a[j][j];
            for (int k = j; k <= dim; ++k) a[j][k] /= d;
            for (int i = 0; i < dim; ++i)
                if (i != j)
                {
                    const double factor = a[i][j];
                    for (int k = j; k <= dim; ++k) a[i][k] -= factor*a[j][k];
                }
        }
        for (int j = 0; j < dim; ++j)
        {
            local[j] += a[j][dim];
            if (!std::isfinite(local[j]) || std::abs(local[j]) > 3) return false;
        }
    }
    return false;
}
}
void RANSTerrainWall::Require(bool valid, const std::string &message) const
{
    int bad = !valid;
    m_session->GetComm()->AllReduce(bad,LibUtilities::ReduceMax);
    if (bad && m_trialMode) throw std::runtime_error("RANS terrain trial: "+message);
    ASSERTL0(!bad,"RANS terrain wall: "+message);
}
NekDouble RANSTerrainWall::Component(const RANSVector &v, int d) const
{
    return v[m_dim == 2 && d == 1 ? 2 : d];
}
RANSVector RANSTerrainWall::Vector(const Fields &v, int q) const
{
    return {v[0][q],m_dim == 3 ? v[1][q] : 0,v[m_dim-1][q]};
}
AWFVector RANSTerrainWall::Tangent(const RANSVector &v, int q) const
{
    return {Dot(v,m_frames[q].tangent),v[1]};
}
RANSVector RANSTerrainWall::Cartesian(const AWFVector &v, int q) const
{
    const auto &t = m_frames[q].tangent;
    return {v[0]*t[0],v[1],v[0]*t[2]};
}

RANSTerrainWall::RANSTerrainWall(const LibUtilities::SessionReaderSharedPtr &session,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields, int dimension,
    const RANSWallGeometry &geometry, const RANSNeutralReference &reference)
    : m_session(session), m_fields(fields), m_dim(dimension), m_terrain(geometry),m_reference(reference)
{
    bool supportedElements = true;
    for (int e = 0; e < fields[0]->GetExpSize(); ++e)
        supportedElements &= fields[0]->GetExp(e)->DetShapeType() ==
            (m_dim == 2 ? LibUtilities::eQuadrilateral : LibUtilities::eHexahedron);
    Require(supportedElements,"terrain wall adapter requires quadrilaterals or hexahedra");
    session->LoadSolverInfo("RANSWallModel",m_wallModel,"None");
    m_log = m_wallModel == "None";
    m_equilibrium = m_wallModel == "AWFEquilibrium";
    m_tag = m_log ? "RANSLogLaw" : m_wallModel;
    m_spatial = ReadRANSWallTransition(session,geometry);
    session->LoadParameter("RANSStripHeight",m_strip.height,2);
    session->LoadParameter("RANSStripPointP",m_strip.pointP,1);
    session->LoadParameter("RANSStripPointN",m_strip.pointN,4);
    m_strip.Validate();
    session->LoadParameter("Kinvis",m_awf.molecularNu);
    std::string resolvedDiffusion;
    session->LoadSolverInfo("RANSResolvedDiffusion",resolvedDiffusion,"Molecular");
    m_resolvedNu = resolvedDiffusion == "Molecular" ? m_awf.molecularNu : 0;
    session->LoadParameter("RANSRoughnessHeight",m_awf.roughness,.1);
    m_awf.height = m_strip.height; m_awf.pointP = m_strip.pointP;
    session->LoadParameter("RANSAWFCmu",m_parameters.cMu,.09);
    session->LoadParameter("RANSAWFCLength",m_parameters.cLength,2.55);
    session->LoadParameter("RANSAWFSmoothIntercept",m_parameters.smoothIntercept,10.7);
    session->LoadParameter("RANSZ0",m_z0,.02);
    session->LoadParameter("RANSHeightShift",m_shift,0);
    session->LoadParameter("RANSMinimumK",m_minimumK,0);
    session->LoadParameter("RANSMinimumEpsilon",m_minimumEpsilon,0);
    session->LoadParameter("RANSNormalAugmentation",m_augmentation,2);
    session->LoadParameter("RANSNormalRelaxation",m_normalRelaxation,.5);
    session->LoadParameter("RANSNormalDirectLimit",m_normalDirectLimit,1024);
    ASSERTL0(m_normalDirectLimit >= 0 && m_normalDirectLimit <= 2048,
             "RANSNormalDirectLimit must be 0..2048 (0 selects iterative projection)");
    session->LoadParameter("RANSMaximumCurvatureRatio",m_maximumCurvature,.25);
    session->LoadParameter("RANSWallGeometryTolerance",m_geometryTolerance,1e-7);
    ASSERTL0(std::isfinite(m_augmentation) && m_augmentation > 0 &&
             std::isfinite(m_normalRelaxation) && m_normalRelaxation > 0 && m_normalRelaxation <= 1 &&
             std::isfinite(m_geometryTolerance) && m_geometryTolerance > 0,
             "Invalid terrain wall constraint parameters");
    m_parameters.Validate();
    RANSWallPropertySignature(session);
    m_failurePolicy = AWFFailurePolicy::UseOriginal;
    if (m_equilibrium)
    {
        ASSERTL0(!session->DefinesParameter("RANSZ0") &&
                 !session->DefinesFunction("RANSWallProperties","z0"),
                 "AWFEquilibrium derives z0 from its matching anchor");
        std::string weight, failure;
        session->LoadSolverInfo("RANSAWFWeight",weight,"Switched");
        session->LoadSolverInfo("RANSAWFFailurePolicy",failure,"UseOriginal");
        ASSERTL0(weight == "Fixed" || weight == "Switched","Invalid RANSAWFWeight");
        ASSERTL0(failure == "UseOriginal" || failure == "RequirePreservation","Invalid AWF failure policy");
        m_switched = weight == "Switched";
        m_failurePolicy = failure == "RequirePreservation" ? AWFFailurePolicy::RequirePreservation : AWFFailurePolicy::UseOriginal;
        session->LoadParameter("RANSAWFFixedWeight",m_fixedWeight,1);
        session->LoadParameter("RANSAWFMinimumRe",m_minimumRe,100);
        ASSERTL0(std::isfinite(m_fixedWeight) && m_fixedWeight >= 0 && m_fixedWeight <= 1 &&
                 std::isfinite(m_minimumRe) && m_minimumRe > 0,"Invalid AWF equilibrium gates");
        if (m_switched)
        {
            session->LoadParameter("RANSAWFOn",m_switch.on);
            session->LoadParameter("RANSAWFOff",m_switch.off);
            const std::array<std::string,7> names{{"U","DeltaU","K","Epsilon","Pressure","Convection","Normal"}};
            for (int i = 0; i < 7; ++i) session->LoadParameter("RANSAWFWeight"+names[i],m_switch.weights[i],1);
            session->LoadParameter("RANSAWFDirectionX",m_switch.incomingDirection[0],1);
            session->LoadParameter("RANSAWFDirectionY",m_switch.incomingDirection[1],0);
            m_switch.Validate();
            ASSERTL0(m_dim != 2 || m_switch.incomingDirection[1] == 0,"2D wall needs streamwise incoming direction");
        }
    }
    for (int b = 0; b < fields[0]->GetBndConditions().size(); ++b)
        if (fields[0]->GetBndConditions()[b]->GetUserDefined() == m_tag)
        {
            m_boundaries.push_back(b); m_offsets.push_back(m_localPoints);
            m_localPoints += fields[0]->GetBndCondExpansions()[b]->GetTotPoints();
        }
    m_x = Array<OneD,NekDouble>(m_localPoints); m_span = Array<OneD,NekDouble>(m_localPoints);
    m_roughness = Array<OneD,NekDouble>(m_localPoints,m_awf.roughness);
    m_logRoughness = Array<OneD,NekDouble>(m_localPoints,m_z0);
    m_mask = Array<OneD,NekDouble>(m_localPoints,1);
    m_frames.resize(m_localPoints); m_points.resize(m_localPoints);
    Fields localCoordinates(m_dim), localShift(m_dim);
    for (auto &c : localCoordinates) c = Array<OneD,NekDouble>(m_localPoints);
    for (auto &c : localShift) c = Array<OneD,NekDouble>(m_localPoints,0.0);
    bool valid = true;
    double geometryDistanceError = 0, geometryNormalError = 0;
    std::string geometryError;
    try
    {
        for (int j = 0; j < m_boundaries.size(); ++j)
        {
            const int b = m_boundaries[j], offset = m_offsets[j];
            auto boundary = fields[0]->GetBndCondExpansions()[b];
            const int count = boundary->GetTotPoints();
            Array<OneD,NekDouble> x(count),y(count),z(count); Fields normals;
            boundary->GetCoords(x,y,z); fields[0]->GetBoundaryNormals(b,normals);
            Array<OneD,NekDouble> px(count),py(count),pz(count),property(count);
            for (int k = 0; k < count; ++k)
            {
                const int q = offset+k;
                const auto f = geometry.Closest(x[k],m_dim == 2 ? y[k] : z[k]);
                geometry.CheckStrip(f,m_strip.height,m_strip.pointN,m_maximumCurvature);
                geometryDistanceError = std::max(geometryDistanceError,std::abs(f.distance-m_strip.height));
                valid &= std::abs(f.distance-m_strip.height) <= m_geometryTolerance*std::max(m_strip.height,.01);
                for (int d = 0; d < m_dim; ++d)
                {
                    geometryNormalError = std::max(geometryNormalError,std::abs(normals[d][k]+Component(f.normal,d)));
                    valid &= std::abs(normals[d][k]+Component(f.normal,d)) <= m_geometryTolerance;
                }
                m_frames[q] = f; m_x[q] = f.x; m_span[q] = m_dim == 2 ? 0 : y[k];
                m_mask[q] = m_spatial.Weight(f.x);
                px[k] = f.x; py[k] = m_dim == 2 ? f.z : y[k]; pz[k] = m_dim == 2 ? 0 : f.z;
                const auto n = geometry.Point(f.x,m_strip.pointN,m_span[q]);
                for (int d = 0; d < m_dim; ++d) localCoordinates[d][q] = Component(n,d);
            }
            // Ownership at a shared tangential vertex belongs to this wall
            // patch's column. Reusing its neighbour's one-sided N derivative
            // would break the integrated strip mass balance.
            for (int e = 0; e < boundary->GetExpSize(); ++e)
            {
                const int start = boundary->GetPhys_Offset(e), size = boundary->GetExp(e)->GetTotPoints();
                double low = m_x[offset+start], high = low, spanLow = m_span[offset+start], spanHigh = spanLow;
                for (int k = start; k < start+size; ++k)
                {
                    low = std::min(low,m_x[offset+k]); high = std::max(high,m_x[offset+k]);
                    spanLow = std::min(spanLow,m_span[offset+k]); spanHigh = std::max(spanHigh,m_span[offset+k]);
                }
                for (int k = start; k < start+size; ++k)
                {
                    const int q = offset+k;
                    if (geometry.terrain == RANSTerrain::CosineHill &&
                        std::abs(std::abs(m_x[q]-geometry.centre)-geometry.hillHalfLength) < 1e-10 &&
                        std::abs(.5*(low+high)-geometry.centre) < geometry.hillHalfLength)
                        m_frames[q].curvature = geometry.hillHeight*std::pow(std::acos(-1.0)/geometry.hillHalfLength,2)/2;
                    for (int d = 0; d < m_dim; ++d)
                        localShift[d][q] = 1e-8*(.5*(low+high)-m_x[q])*Component(m_frames[q].tangent,d);
                    if (m_dim == 3) localShift[1][q] = 1e-8*(.5*(spanLow+spanHigh)-m_span[q]);
                }
            }
            if (session->DefinesFunction("RANSWallProperties","roughness"))
            {
                session->GetFunction("RANSWallProperties","roughness")->Evaluate(px,py,pz,0,property);
                Vmath::Vcopy(count,property.data(),1,m_roughness.data()+offset,1);
            }
            if (session->DefinesFunction("RANSWallProperties","z0"))
            {
                session->GetFunction("RANSWallProperties","z0")->Evaluate(px,py,pz,0,property);
                Vmath::Vcopy(count,property.data(),1,m_logRoughness.data()+offset,1);
            }
        }
        for (int q = 0; q < m_localPoints; ++q)
        {
            auto g = m_awf; g.roughness = m_roughness[q];
            if (!m_log) g.Validate();
            valid &= !m_log || (std::isfinite(m_logRoughness[q]) && m_logRoughness[q] > 0 &&
                                m_strip.height+m_shift > m_logRoughness[q]);
            valid &= !m_equilibrium || m_mask[q] == 0 || (!m_switched && m_fixedWeight == 0) ||
                     std::abs(m_frames[q].curvature) < 1e-12;
        }
    }
    catch (const std::exception &error) { valid = false; geometryError = error.what(); }
    if (!valid)
        std::cerr << std::setprecision(17) << "Terrain geometry rank " << session->GetComm()->GetRank()
                  << " distance error " << geometryDistanceError << " normal error " << geometryNormalError
                  << " detail: " << geometryError << std::endl;
    Require(valid,"interface must be a physical normal offset of the terrain; invalid wall properties/curvature, or equilibrium correction active on curved terrain");
    auto comm = session->GetComm();
    Array<OneD,int> count(1,m_localPoints), counts, offsets(comm->GetSize());
    comm->AllGather(count,counts);
    for (int r = 0; r < counts.size(); ++r)
    {
        offsets[r] = m_globalPoints;
        if (r == comm->GetRank()) m_globalOffset = m_globalPoints;
        m_globalPoints += counts[r];
    }
    Require(m_globalPoints > 0,"no matching terrain wall regions");
    Fields coordinates(m_dim), shifts(m_dim);
    for (int d = 0; d < m_dim; ++d)
    {
        coordinates[d] = Array<OneD,NekDouble>(m_globalPoints);
        comm->AllGatherv(localCoordinates[d],coordinates[d],counts,offsets);
        shifts[d] = Array<OneD,NekDouble>(m_globalPoints);
        comm->AllGatherv(localShift[d],shifts[d],counts,offsets);
    }
    m_sampleElements.resize(m_globalPoints,-1); m_sampleCoordinates.resize(m_globalPoints);
    Array<OneD,int> ids(m_globalPoints,std::numeric_limits<int>::max());
    for (int q = 0; q < m_globalPoints; ++q)
    {
        Array<OneD,NekDouble> point(m_dim), local(m_dim);
        for (int d = 0; d < m_dim; ++d) point[d] = coordinates[d][q];
        const auto f = geometry.Closest(point[0],point[m_dim-1]);
        // Move only ownership slightly above N; evaluate the exact physical N.
        for (int d = 0; d < m_dim; ++d) point[d] += shifts[d][q]+1e-8*m_strip.height*Component(f.normal,d);
        int e = fields[0]->GetExpIndex(point,local,1e-7);
        double residual = 0;
        if (e >= 0 && !RefineCoordinates(fields[0]->GetExp(e)->GetGeom(),point,local,residual,1e-10)) e = -1;
        if (e < 0)
            for (int candidate = 0; candidate < fields[0]->GetExpSize(); ++candidate)
            {
                const auto geom = fields[0]->GetExp(candidate)->GetGeom();
                const auto box = geom->GetBoundingBox(); bool inBox = true;
                for (int d = 0; d < m_dim; ++d) inBox &= point[d] >= box[d]-1e-8 && point[d] <= box[d+3]+1e-8;
                if (inBox && RefineCoordinates(geom,point,local,residual,1e-10)) { e = candidate; break; }
            }
        if (e >= 0)
        {
            for (int d = 0; d < m_dim; ++d) point[d] = coordinates[d][q];
            if (!RefineCoordinates(fields[0]->GetExp(e)->GetGeom(),point,local,residual,1e-7)) continue;
            m_samplingResidual = std::max(m_samplingResidual,residual);
            m_sampleElements[q] = e; m_sampleCoordinates[q] = local;
            ids[q] = fields[0]->GetExp(e)->GetGeom()->GetGlobalID();
        }
    }
    comm->AllReduce(ids,LibUtilities::ReduceMin);
    Array<OneD,int> owners(m_globalPoints,std::numeric_limits<int>::max());
    valid = true;
    for (int q = 0; q < m_globalPoints; ++q)
    {
        valid &= ids[q] != std::numeric_limits<int>::max();
        const int e = m_sampleElements[q];
        if (e >= 0 && fields[0]->GetExp(e)->GetGeom()->GetGlobalID() == ids[q]) owners[q] = comm->GetRank();
        else m_sampleElements[q] = -1;
    }
    Require(valid,"every normal N sample must be inside the resolved domain");
    comm->AllReduce(m_samplingResidual,LibUtilities::ReduceMax);
    comm->AllReduce(owners,LibUtilities::ReduceMin);
    for (int q = 0; q < m_globalPoints; ++q)
        if (m_sampleElements[q] >= 0 && owners[q] != comm->GetRank()) m_sampleElements[q] = -1;
    for (int q = 0; q < m_localPoints; ++q) m_remoteSamples += owners[m_globalOffset+q] != comm->GetRank();
    comm->AllReduce(m_remoteSamples,LibUtilities::ReduceSum);
    for (auto *a : {&m_multiplier,&m_normalTarget,&m_normalError,&m_projectedError,&m_velocityMass})
        *a = Array<OneD,NekDouble>(m_localPoints,0.0);
    m_diagonal = Fields(m_dim+1); m_rhs = Fields(m_dim+1);
    for (int i = 0; i <= m_dim; ++i)
    {
        m_diagonal[i] = Array<OneD,NekDouble>(m_localPoints,0.0);
        m_rhs[i] = Array<OneD,NekDouble>(m_localPoints,0.0);
    }
    // A continuous velocity trace avoids both the uncontrolled high modes of
    // the pressure trace and redundant face-endpoint constraints of a DG trace.
    // Tests at prescribed-velocity corners are excluded, as in the volume CG
    // equations. Assemble the trace mass diagonal using the native CG map.
    const auto map = std::dynamic_pointer_cast<MultiRegions::ContField>(m_fields[0])->GetLocalToGlobalMap();
    m_traceOffsets.resize(m_fields[0]->GetBndCondExpansions().size());
    int coefficientOffset = 0;
    for (int b = 0; b < m_traceOffsets.size(); ++b)
    {
        m_traceOffsets[b] = coefficientOffset;
        coefficientOffset += m_fields[0]->GetBndCondExpansions()[b]->GetNcoeffs();
    }
    Array<OneD,NekDouble> diagonal(m_fields[0]->GetNcoeffs(),0.0);
    for (int b : m_boundaries)
    {
        auto boundary = m_fields[0]->GetBndCondExpansions()[b];
        for (int e = 0; e < boundary->GetExpSize(); ++e)
        {
            auto expansion = boundary->GetExp(e);
            const auto mass = expansion->GetLocMatrix(StdRegions::eMass);
            const int start = m_traceOffsets[b]+boundary->GetCoeff_Offset(e);
            for (int k = 0; k < expansion->GetNcoeffs(); ++k)
            {
                const int local = map->GetBndCondCoeffsToLocalCoeffsMap()[start+k];
                diagonal[local] += map->GetLocalToGlobalSign(local)*(*mass)(k,k);
            }
        }
    }
    m_normalDiagonal = Array<OneD,NekDouble>(map->GetNumGlobalCoeffs());
    map->Assemble(diagonal,m_normalDiagonal);
    for (int i = 0; i < map->GetNumGlobalDirBndCoeffs(); ++i) m_normalDiagonal[i] = 0;
    BuildNormalDirect();
    Array<OneD,NekDouble> a(m_localPoints),b(m_localPoints),work(m_localPoints),ones(m_localPoints,1.0);
    for (int q = 0; q < m_localPoints; ++q)
    {
        a[q] = std::sin(37*m_x[q]+19*m_span[q]);
        b[q] = std::cos(23*m_x[q]-17*m_span[q]);
    }
    const auto pa = ProjectNormal(a), pb = ProjectNormal(b), ppa = ProjectNormal(pa);
    for (int q = 0; q < m_localPoints; ++q) work[q] = a[q]*pb[q]-pa[q]*b[q];
    const double area = SurfaceIntegral(ones);
    m_projectionAudit = std::abs(SurfaceIntegral(work))/std::max(area,1e-30);
    for (int q = 0; q < m_localPoints; ++q) work[q] = (pa[q]-ppa[q])*(pa[q]-ppa[q]);
    m_projectionAudit = std::max(m_projectionAudit,std::sqrt(std::max(0.0,SurfaceIntegral(work))/area));
    const auto constant = ProjectNormal(ones);
    for (int q = 0; q < m_localPoints; ++q) work[q] = (a[q]-pa[q])*constant[q];
    m_projectionAudit = std::max(m_projectionAudit,std::abs(SurfaceIntegral(work))/area);
    m_constantProjection = constant;
    m_constantComplement = Array<OneD,NekDouble>(m_localPoints);
    for (int q = 0; q < m_localPoints; ++q) m_constantComplement[q] = 1-constant[q];
    m_constraintMeasure = SurfaceIntegral(constant);
    Require(m_constraintMeasure > 0,"empty free normal trace space");
    const auto ca = ProjectConstraint(a),cb = ProjectConstraint(b),cca = ProjectConstraint(ca);
    const auto adjoint = AdjointConstraint(pb);
    for (int q = 0; q < m_localPoints; ++q) work[q] = ca[q]*pb[q]-a[q]*adjoint[q];
    m_projectionAudit = std::max(m_projectionAudit,std::abs(SurfaceIntegral(work))/area);
    for (int q = 0; q < m_localPoints; ++q) work[q] = (ca[q]-cca[q])*(ca[q]-cca[q]);
    m_projectionAudit = std::max(m_projectionAudit,std::sqrt(std::max(0.0,SurfaceIntegral(work))/area));
    for (int q = 0; q < m_localPoints; ++q) work[q] = ca[q]-a[q];
    m_projectionAudit = std::max(m_projectionAudit,std::abs(SurfaceIntegral(work))/area);
    const auto aca = AdjointConstraint(ca),acb = AdjointConstraint(cb);
    for (int q = 0; q < m_localPoints; ++q) work[q] = a[q]*acb[q]-b[q]*aca[q];
    m_projectionAudit = std::max(m_projectionAudit,std::abs(SurfaceIntegral(work))/area);
    Require(m_projectionAudit < 5e-10,"invalid normal projector / conservative constraint adjoint");
    for (int q = 0; q < m_localPoints; ++q)
        work[q] = m_frames[q].Volume(m_strip.height)/m_frames[q].Jacobian(m_strip.height);
    m_unresolvedVolume = SurfaceIntegral(work);
    Array<OneD,NekDouble> volumeOne(m_fields[0]->GetTotPoints(),1.0);
    m_resolvedVolume = m_fields[0]->Integral(volumeOne);
    if (comm->GetRank() == 0)
    {
        m_diagnostics.open(session->GetSessionName()+".terrain-wall.csv");
        m_diagnostics << "time,normal_projected_max,normal_pointwise_max,epsilon_projection_max,wall_points,remote_samples,chi_min,chi_max,min_kP,min_epsH,area,physical_wall_area,N_inverse_residual,projection_audit,resolved_volume,unresolved_volume,normal_integral_defect,resolved_mass_defect,tangent_flux_residual_l2,k_flux_residual_l2,normal_l2,normal_away_feet_max\n" << std::setprecision(17);
    }
    int write = 0; session->LoadParameter("RANSAWFPointDiagnostics",write,0);
    ASSERTL0(write == 0 || write == 1,"RANSAWFPointDiagnostics must be 0 or 1");
    if (write)
    {
        m_pointDiagnostics.open(session->GetSessionName()+".terrain-points-r"+std::to_string(comm->GetRank())+".csv");
        m_pointDiagnostics << "time,phase,point,x,span,curvature,JH,V,H,P,N,roughness,z0,mask,kP,uH,uP,uN,kN,epsilonN,forcing_x,forcing_span,traction_x,traction_span,production,dissipation,epsP,epsH,chi,requested_chi,departure,reason,wall_re,cutoff,normal_target,normal_error,normal_projected_error,multiplier,tangent_source,k_source,spanH,spanP,spanN,pressure_t,pressure_span,convection_t,convection_span,wP,r_U,r_DeltaU,r_K,r_Epsilon,r_Pressure,r_Convection,r_Normal,sensors_evaluated,resolved_shear,tangent_flux_residual,k_flux_residual\n" << std::setprecision(17);
    }
}

RANSTerrainWall::Fields RANSTerrainWall::Trace(const Fields &v) const
{
    Fields result(v.size());
    for (int f = 0; f < v.size(); ++f)
    {
        result[f] = Array<OneD,NekDouble>(m_localPoints);
        for (int j = 0; j < m_boundaries.size(); ++j)
        {
            Array<OneD,NekDouble> trace;
            m_fields[f]->ExtractPhysToBnd(m_boundaries[j],v[f],trace);
            Vmath::Vcopy(trace.size(),trace.data(),1,result[f].data()+m_offsets[j],1);
        }
    }
    return result;
}
RANSTerrainWall::Fields RANSTerrainWall::SampleN(const Fields &v) const
{
    Fields result(v.size());
    for (int f = 0; f < v.size(); ++f)
    {
        Array<OneD,NekDouble> all(m_globalPoints,0.0);
        for (int q = 0; q < m_globalPoints; ++q)
        {
            const int e = m_sampleElements[q];
            if (e >= 0) all[q] = m_fields[0]->GetExp(e)->StdPhysEvaluate(
                m_sampleCoordinates[q],v[f]+m_fields[0]->GetPhys_Offset(e));
        }
        m_session->GetComm()->AllReduce(all,LibUtilities::ReduceSum);
        result[f] = Array<OneD,NekDouble>(m_localPoints);
        Vmath::Vcopy(m_localPoints,all+m_globalOffset,1,result[f],1);
    }
    return result;
}
void RANSTerrainWall::SetBoundary(int field, const Array<OneD,NekDouble> &values, bool dirichlet)
{
    for (int j = 0; j < m_boundaries.size(); ++j)
    {
        auto b = m_fields[field]->GetBndCondExpansions()[m_boundaries[j]];
        Array<OneD,NekDouble> local = values+m_offsets[j];
        if (dirichlet) b->FwdTransBndConstrained(local,b->UpdateCoeffs());
        else b->IProductWRTBase(local,b->UpdateCoeffs());
    }
}
void RANSTerrainWall::SetIterationState(const Array<OneD,NekDouble> &state)
{
    ASSERTL0(state.size() == m_multiplier.size(),"terrain iteration state size mismatch");
    m_multiplier = ProjectNormal(state);
}
Array<OneD,NekDouble> RANSTerrainWall::IterationMetric() const
{
    Array<OneD,NekDouble> result(m_localPoints);
    for (int j = 0; j < m_boundaries.size(); ++j)
    {
        auto boundary = m_fields[0]->GetBndCondExpansions()[m_boundaries[j]];
        for (int e = 0; e < boundary->GetExpSize(); ++e)
        {
            auto expansion = boundary->GetExp(e);
            const int offset = m_offsets[j]+boundary->GetPhys_Offset(e);
            Array<OneD,NekDouble> one(expansion->GetTotPoints(),1.0),out = result+offset;
            expansion->MultiplyByQuadratureMetric(one,out);
        }
    }
    for (int q = 0; q < m_localPoints; ++q)
        result[q] *= m_frames[q].Volume(m_strip.height)/m_frames[q].Jacobian(m_strip.height);
    return result;
}
Array<OneD,NekDouble> RANSTerrainWall::ProjectConstraint(const Array<OneD,NekDouble> &values) const
{
    auto result = ProjectNormal(values);
    Array<OneD,NekDouble> difference(values.size());
    for (int q = 0; q < values.size(); ++q) difference[q] = values[q]-result[q];
    const double correction = SurfaceIntegral(difference)/m_constraintMeasure;
    for (int q = 0; q < values.size(); ++q) result[q] += correction*m_constantProjection[q];
    return result;
}
Array<OneD,NekDouble> RANSTerrainWall::AdjointConstraint(const Array<OneD,NekDouble> &values) const
{
    // Input belongs to the free continuous trace. C* lifts its dual tests
    // toward the constant, restoring global strip conservation at prescribed
    // inlet/wall corners. C* C is self-adjoint positive semidefinite.
    Array<OneD,NekDouble> result(values.size()),product(values.size());
    for (int q = 0; q < values.size(); ++q) product[q] = m_constantProjection[q]*values[q];
    const double correction = SurfaceIntegral(product)/m_constraintMeasure;
    for (int q = 0; q < values.size(); ++q) result[q] = values[q]+correction*m_constantComplement[q];
    return result;
}
Array<OneD,NekDouble> RANSTerrainWall::ProjectNormal(const Array<OneD,NekDouble> &values) const
{
    const auto map = std::dynamic_pointer_cast<MultiRegions::ContField>(m_fields[0])->GetLocalToGlobalMap();
    auto load = NormalLoad(values);
    if (m_normalDirectSize)
    {
        const int n = m_normalDirectSize;
        Array<OneD,NekDouble> rhs(n,0.0),solution(load.size(),0.0);
        for (int i = 0; i < load.size(); ++i)
            if (m_normalDirectMap[i] >= 0 && map->GetGlobalToUniversalMapUnique(i) > 0)
                rhs[m_normalDirectMap[i]] = load[i];
        m_session->GetComm()->AllReduce(rhs,LibUtilities::ReduceSum);
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < i; ++j) rhs[i] -= m_normalCholesky[i*n+j]*rhs[j];
            rhs[i] /= m_normalCholesky[i*n+i];
        }
        for (int i = n-1; i >= 0; --i)
        {
            for (int j = i+1; j < n; ++j) rhs[i] -= m_normalCholesky[j*n+i]*rhs[j];
            rhs[i] /= m_normalCholesky[i*n+i];
        }
        for (int i = 0; i < solution.size(); ++i)
            if (m_normalDirectMap[i] >= 0) solution[i] = rhs[m_normalDirectMap[i]];
        return NormalValues(solution);
    }
    const int count = load.size();
    Array<OneD,NekDouble> solution(count,0.0),residual(count),z(count),direction(count),applied(count);
    Vmath::Vcopy(count,load,1,residual,1);
    auto dot = [&](const Array<OneD,NekDouble> &a, const Array<OneD,NekDouble> &b) {
        double sum = 0;
        for (int i = 0; i < count; ++i)
            if (m_normalDiagonal[i] > 0 && map->GetGlobalToUniversalMapUnique(i) > 0) sum += a[i]*b[i];
        m_session->GetComm()->AllReduce(sum,LibUtilities::ReduceSum);
        return sum;
    };
    const double norm = std::sqrt(std::max(0.0,dot(load,load)));
    const double tolerance = std::max(2e-16,2e-15*norm);
    double previous = 0;
    int iteration = 0;
    for (; iteration < 300; ++iteration)
    {
        if (std::sqrt(std::max(0.0,dot(residual,residual))) <= tolerance) break;
        for (int i = 0; i < count; ++i) z[i] = m_normalDiagonal[i] > 0 ?
            residual[i]/m_normalDiagonal[i] : 0;
        const double next = dot(residual,z);
        ASSERTL0(std::isfinite(next) && next > 0,"terrain CG trace mass is not positive definite");
        const double beta = iteration ? next/previous : 0;
        for (int i = 0; i < count; ++i) direction[i] = z[i]+beta*(iteration ? direction[i] : 0);
        applied = NormalLoad(NormalValues(direction));
        const double curvature = dot(direction,applied);
        ASSERTL0(std::isfinite(curvature) && curvature > 0,"terrain CG trace mass curvature is not positive");
        const double step = next/curvature;
        for (int i = 0; i < count; ++i)
        {
            solution[i] += step*direction[i]; residual[i] -= step*applied[i];
        }
        previous = next;
    }
    ASSERTL0(iteration < 300,"terrain CG trace projection failed to converge");
    return NormalValues(solution);
}
Array<OneD,NekDouble> RANSTerrainWall::NormalLoad(const Array<OneD,NekDouble> &values) const
{
    const auto map = std::dynamic_pointer_cast<MultiRegions::ContField>(m_fields[0])->GetLocalToGlobalMap();
    Array<OneD,NekDouble> local(m_fields[0]->GetNcoeffs(),0.0),global(map->GetNumGlobalCoeffs());
    for (int j = 0; j < m_boundaries.size(); ++j)
    {
        auto b = m_fields[0]->GetBndCondExpansions()[m_boundaries[j]];
        Array<OneD,NekDouble> coefficients(b->GetNcoeffs());
        b->IProductWRTBase(values+m_offsets[j],coefficients);
        const int start = m_traceOffsets[m_boundaries[j]];
        for (int i = 0; i < coefficients.size(); ++i)
            local[map->GetBndCondCoeffsToLocalCoeffsMap()[start+i]] +=
                map->GetBndCondCoeffsToLocalCoeffsSign()[start+i]*coefficients[i];
    }
    map->Assemble(local,global);
    for (int i = 0; i < global.size(); ++i) if (m_normalDiagonal[i] <= 0) global[i] = 0;
    return global;
}
Array<OneD,NekDouble> RANSTerrainWall::NormalValues(const Array<OneD,NekDouble> &global) const
{
    const auto map = std::dynamic_pointer_cast<MultiRegions::ContField>(m_fields[0])->GetLocalToGlobalMap();
    Array<OneD,NekDouble> local(m_fields[0]->GetNcoeffs()),result(m_localPoints);
    map->GlobalToLocal(global,local);
    for (int j = 0; j < m_boundaries.size(); ++j)
    {
        auto b = m_fields[0]->GetBndCondExpansions()[m_boundaries[j]];
        const int start = m_traceOffsets[m_boundaries[j]];
        Array<OneD,NekDouble> coefficients(b->GetNcoeffs()),out = result+m_offsets[j];
        for (int i = 0; i < coefficients.size(); ++i)
            coefficients[i] = map->GetBndCondCoeffsToLocalCoeffsSign()[start+i]*
                local[map->GetBndCondCoeffsToLocalCoeffsMap()[start+i]];
        b->BwdTrans(coefficients,out);
    }
    return result;
}
void RANSTerrainWall::BuildNormalDirect()
{
    if (!m_normalDirectLimit) return;
    const auto map = std::dynamic_pointer_cast<MultiRegions::ContField>(m_fields[0])->GetLocalToGlobalMap();
    auto comm = m_session->GetComm();
    std::vector<int> localIds;
    for (int i = 0; i < m_normalDiagonal.size(); ++i)
        if (m_normalDiagonal[i] > 0 && map->GetGlobalToUniversalMapUnique(i) > 0)
            localIds.push_back(map->GetGlobalToUniversalMap()[i]);
    Array<OneD,int> count(1,localIds.size()),counts,offsets(comm->GetSize());
    comm->AllGather(count,counts);
    int size = 0;
    for (int r = 0; r < counts.size(); ++r) { offsets[r] = size; size += counts[r]; }
    if (!size || size > m_normalDirectLimit) return;
    Array<OneD,int> ids(localIds.size()),all(size);
    for (int i = 0; i < ids.size(); ++i) ids[i] = localIds[i];
    comm->AllGatherv(ids,all,counts,offsets);
    std::vector<int> ordered(all.begin(),all.end());
    std::sort(ordered.begin(),ordered.end());
    ordered.erase(std::unique(ordered.begin(),ordered.end()),ordered.end());
    const int n = ordered.size();
    m_normalDirectMap.resize(m_normalDiagonal.size(),-1);
    for (int i = 0; i < m_normalDiagonal.size(); ++i)
        if (m_normalDiagonal[i] > 0)
        {
            auto found = std::lower_bound(ordered.begin(),ordered.end(),map->GetGlobalToUniversalMap()[i]);
            ASSERTL0(found != ordered.end() && *found == map->GetGlobalToUniversalMap()[i],"missing global terrain trace mode");
            m_normalDirectMap[i] = found-ordered.begin();
        }
    Array<OneD,NekDouble> matrix(n*n,0.0);
    for (int b : m_boundaries)
    {
        auto boundary = m_fields[0]->GetBndCondExpansions()[b];
        for (int e = 0; e < boundary->GetExpSize(); ++e)
        {
            auto expansion = boundary->GetExp(e);
            const auto mass = expansion->GetLocMatrix(StdRegions::eMass);
            const int start = m_traceOffsets[b]+boundary->GetCoeff_Offset(e);
            for (int i = 0; i < expansion->GetNcoeffs(); ++i)
                for (int j = 0; j < expansion->GetNcoeffs(); ++j)
                {
                    const int li = map->GetBndCondCoeffsToLocalCoeffsMap()[start+i],lj = map->GetBndCondCoeffsToLocalCoeffsMap()[start+j];
                    const int gi = m_normalDirectMap[map->GetLocalToGlobalMap(li)],gj = m_normalDirectMap[map->GetLocalToGlobalMap(lj)];
                    if (gi < 0 || gj < 0) continue;
                    const double si = map->GetBndCondCoeffsToLocalCoeffsSign()[start+i]*map->GetLocalToGlobalSign(li);
                    const double sj = map->GetBndCondCoeffsToLocalCoeffsSign()[start+j]*map->GetLocalToGlobalSign(lj);
                    matrix[gi*n+gj] += si*sj*(*mass)(i,j);
                }
        }
    }
    comm->AllReduce(matrix,LibUtilities::ReduceSum);
    m_normalCholesky.assign(n*n,0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j <= i; ++j)
        {
            double value = matrix[i*n+j];
            for (int k = 0; k < j; ++k) value -= m_normalCholesky[i*n+k]*m_normalCholesky[j*n+k];
            ASSERTL0(std::isfinite(value) && (i != j || value > 0),"terrain trace mass Cholesky factorization failed");
            m_normalCholesky[i*n+j] = i == j ? std::sqrt(value) : value/m_normalCholesky[j*n+j];
        }
    m_normalDirectSize = n;
}
NekDouble RANSTerrainWall::SurfaceIntegral(const Array<OneD,NekDouble> &values) const
{
    double result = 0;
    for (int j = 0; j < m_boundaries.size(); ++j)
    {
        auto boundary = m_fields[0]->GetBndCondExpansions()[m_boundaries[j]];
        for (int e = 0; e < boundary->GetExpSize(); ++e)
            result += boundary->GetExp(e)->Integral(values+m_offsets[j]+boundary->GetPhys_Offset(e));
    }
    m_session->GetComm()->AllReduce(result,LibUtilities::ReduceSum);
    return result;
}

void RANSTerrainWall::CheckInitialState(const Fields &state) const
{
    if (!m_equilibrium || m_failurePolicy != AWFFailurePolicy::RequirePreservation) return;
    const auto h = Trace(state), n = SampleN(state);
    bool valid = true;
    try
    {
        for (int q = 0; q < m_localPoints; ++q)
        {
            if (m_mask[q] == 0 || (!m_switched && m_fixedWeight == 0)) continue;
            const double k = m_strip.EliminateKP(h[m_dim][q],n[m_dim][q]);
            auto g = m_awf; g.roughness = m_roughness[q];
            const EquilibriumWallFunction wall(g,k,m_reference,m_parameters);
            valid &= wall.Reference().admissible && std::sqrt(k)*g.height/g.molecularNu >= m_minimumRe;
        }
    }
    catch (const std::exception &) { valid = false; }
    Require(valid,"requested equilibrium reference fails cutoff/high-Re preflight");
}

std::string RANSTerrainWall::RestartSignature() const
{
    std::ostringstream s; s << std::setprecision(17) << "TerrainWall-10:conservative-dual-surface-continuity:" << m_wallModel << ':'
        << static_cast<int>(m_terrain.terrain) << ':' << m_terrain.origin << ':' << m_terrain.centre << ':'
        << m_terrain.slope << ':' << m_terrain.hillHeight << ':' << m_terrain.hillHalfLength << ':'
        << m_strip.height << ':' << m_strip.pointP << ':' << m_strip.pointN << ':' << m_awf.roughness << ':'
        << m_parameters.cMu << ':' << m_parameters.cLength << ':' << m_parameters.smoothIntercept << ':'
        << m_spatial.enabled << ':' << m_spatial.begin << ':' << m_spatial.end << ':' << m_switched << ':'
        << m_fixedWeight << ':' << m_minimumRe << ':' << static_cast<int>(m_failurePolicy) << ':'
        << m_augmentation << ':' << m_normalRelaxation << ':' << m_normalDirectLimit << ':' << m_maximumCurvature << ':' << m_geometryTolerance;
    if (m_switched)
    {
        s << ':' << m_switch.on << ':' << m_switch.off;
        for (double w : m_switch.weights) s << ':' << w;
        for (double d : m_switch.incomingDirection) s << ':' << d;
    }
    s << RANSWallPropertySignature(m_session);
    return s.str();
}

void RANSTerrainWall::UpdateContinuity(const Fields &h,const Fields &n)
{
    m_streamDivergence = Array<OneD,NekDouble>(m_localPoints,0.0);
    m_spanDivergence = Array<OneD,NekDouble>(m_localPoints,0.0);
    if (m_log) return;
    // Differentiate the pulled-back scalar strip flux, after rotating H/N
    // velocities into the physical wall frame. This is conservative on the
    // discrete interface. Differentiating volume gradients and then applying
    // a continuous chain rule aliases the rotating frame on a curved mesh.
    for (int j = 0; j < m_boundaries.size(); ++j)
    {
        auto boundary = m_fields[0]->GetBndCondExpansions()[m_boundaries[j]];
        const int count = boundary->GetTotPoints(),offset = m_offsets[j];
        Array<OneD,NekDouble> stream(count),span(count),derivative(count);
        for (int k = 0; k < count; ++k)
        {
            const int q = offset+k;
            const auto uh = Tangent(Vector(h,q),q),un = Tangent(Vector(n,q),q);
            stream[k] = m_strip.EliminateVelocityP(uh[0],un[0]);
            span[k] = m_strip.EliminateVelocityP(uh[1],un[1]);
        }
        for (int i = 0; i < m_dim; ++i)
        {
            boundary->PhysDeriv(MultiRegions::DirCartesianMap[i],stream,derivative);
            for (int k = 0; k < count; ++k)
                m_streamDivergence[offset+k] += m_frames[offset+k].Jacobian(m_strip.height)*
                    Component(m_frames[offset+k].tangent,i)*derivative[k];
        }
        if (m_dim == 3)
        {
            boundary->PhysDeriv(MultiRegions::DirCartesianMap[1],span,derivative);
            Array<OneD,NekDouble> out = m_spanDivergence+offset;
            Vmath::Vcopy(count,derivative,1,out,1);
        }
    }
}

void RANSTerrainWall::Update(const Fields &state, const Fields &predictor, const Gradients &gradient,
    const Fields &pressureForce, const Fields &convection, const Fields &bodyForce, NekDouble alpha)
{
    const auto h = Trace(state), n = SampleN(state), bh = Trace(predictor), bn = SampleN(predictor);
    const auto ph = Trace(pressureForce), pn = SampleN(pressureForce);
    const auto ah = Trace(convection), an = SampleN(convection), fh = Trace(bodyForce), fn = SampleN(bodyForce);
    Fields derivatives(m_dim*m_dim);
    for (int i = 0; i < m_dim; ++i) for (int j = 0; j < m_dim; ++j) derivatives[i*m_dim+j] = gradient[i][j];
    // Derivatives belong to velocity fields; Trace uses the corresponding
    // matching expansions, whose ordered boundary grid is shared by all fields.
    auto sampleDerivatives = [&](bool remote) {
        Fields result(m_dim*m_dim);
        for (int i = 0; i < result.size(); ++i)
        {
            Fields one(1); one[0] = derivatives[i];
            result[i] = (remote ? SampleN(one) : Trace(one))[0];
        }
        return result;
    };
    const auto dh = sampleDerivatives(false);
    UpdateContinuity(h,n);
    if (!m_multiplierInitialized)
    {
        Array<OneD,NekDouble> reaction(m_localPoints,0.0);
        bool admissible = true;
        try
        {
            for (int q = 0; q < m_localPoints; ++q)
            {
                RANSTensor g{};
                for (int i = 0; i < m_dim; ++i)
                    for (int j = 0; j < m_dim; ++j)
                        g[m_dim == 2 && i == 1 ? 2 : i][m_dim == 2 && j == 1 ? 2 : j] = dh[i*m_dim+j][q];
                const auto c = EvaluateRANSClosure(m_reference.Model(),h[m_dim][q],h[m_dim+1][q],
                    g,m_strip.height,m_reference.Parameters(),m_resolvedNu);
                const auto &n = m_frames[q].normal;
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < 3; ++j)
                        reaction[q] -= n[i]*(c.momentumDiffusion*(g[i][j]+g[j][i])-
                            c.quadraticStress[i][j]-c.cubicStress[i][j])*n[j];
            }
        }
        catch (const std::exception &) { admissible = false; }
        Require(admissible,"cannot initialize the wall-normal stress reaction");
        m_multiplier = ProjectNormal(reaction);
        m_multiplierInitialized = true;
    }
    const double w = m_strip.VelocityWeightN(), wk = m_strip.KWeightN();
    Array<OneD,NekDouble> beta(m_localPoints), eps(m_localPoints);
    std::vector<AWFVector> tangentRhs(m_localPoints);
    Array<OneD,NekDouble> tangentMass(m_localPoints), kMass(m_localPoints), kRhs(m_localPoints);
    bool valid = true;
    try
    {
        for (int q = 0; q < m_localPoints; ++q)
        {
            auto &p = m_points[q]; p = Point{};
            const auto &frame = m_frames[q]; const auto &normal = frame.normal;
            const double J = frame.Jacobian(m_strip.height), V = frame.Volume(m_strip.height);
            const auto uh = Vector(h,q), un = Vector(n,q);
            const auto th = Tangent(uh,q), tn = Tangent(un,q);
            const auto bth = Tangent(Vector(bh,q),q), btn = Tangent(Vector(bn,q),q);
            // Integrated strip continuity per physical wall area:
            // J(H) W_H + H d_s(U_P) + V(H) d_span(U_span,P) = 0.
            m_normalTarget[q] = m_log ? 0 :
                -(m_strip.height*m_streamDivergence[q]+V*m_spanDivergence[q])/J;
            m_normalError[q] = Dot(uh,normal)-m_normalTarget[q];
            p.kP = m_log ? h[m_dim][q] : m_strip.EliminateKP(h[m_dim][q],n[m_dim][q]);
            if (!(p.kP > m_minimumK)) throw std::invalid_argument("nonpositive terrain KP");
            beta[q] = m_augmentation*(alpha*V/J+std::sqrt(p.kP));
            if (m_log)
            {
                const auto tau = m_reference.StandardLogWallTraction(p.kP,th,m_strip.height,m_logRoughness[q],m_shift);
                const double us = std::sqrt(p.kP/m_reference.Q());
                eps[q] = us*us*us/(m_reference.Kappa()*(m_strip.height+m_shift));
                tangentMass[q] = 0; tangentRhs[q] = {-tau[0],-tau[1]}; kMass[q] = 0; kRhs[q] = 0;
                p.selected.values.traction = tau; p.selected.values.dissipationH = eps[q];
                p.epsH = eps[q];
                continue;
            }
            auto g = m_awf; g.roughness = m_roughness[q];
            const AnalyticalWallFunction wall(g,p.kP,m_parameters);
            auto &sample = p.sample;
            sample.matchingVelocity = th; sample.heightN = m_strip.pointN;
            sample.kN = n[m_dim][q]; sample.epsilonN = n[m_dim+1][q];
            sample.normalVelocityP = -(m_strip.pointP*m_streamDivergence[q]+
                frame.Volume(m_strip.pointP)*m_spanDivergence[q])/frame.Jacobian(m_strip.pointP);
            const auto pressureH = Tangent(Vector(ph,q),q), pressureN = Tangent(Vector(pn,q),q);
            const auto bodyH = Tangent(Vector(fh,q),q), bodyN = Tangent(Vector(fn,q),q);
            const auto advH = Tangent(Vector(ah,q),q), advN = Tangent(Vector(an,q),q);
            AWFVector force{}, pressure{};
            for (int i = 0; i < 2; ++i)
            {
                sample.velocityP[i] = m_strip.EliminateVelocityP(th[i],tn[i]); sample.velocityN[i] = tn[i];
                pressure[i] = m_strip.EliminateVelocityP(pressureH[i],pressureN[i]);
                sample.pressureForcing[i] = pressure[i]-m_strip.EliminateVelocityP(bodyH[i],bodyN[i]);
                sample.convectionForcing[i] = m_strip.EliminateVelocityP(advH[i],advN[i]);
                force[i] = sample.pressureForcing[i]+sample.convectionForcing[i];
            }
            p.selected.values = wall.Evaluate(th,force);
            p.selected.reason = AWFFallbackReason::SwitchInactive;
            p.wallRe = std::sqrt(p.kP)*g.height/g.molecularNu;
            double numerator = 0, cutoff = 0;
            if (m_equilibrium && m_mask[q] > 0 && (m_switched || m_fixedWeight > 0))
            {
                const EquilibriumWallFunction corrected(g,p.kP,m_reference,m_parameters);
                const auto &r = corrected.Reference(); numerator = r.numerator; cutoff = r.cutoff; p.cutoff = cutoff;
                if (m_switched)
                {
                    // Obtain sensors first, then apply the spatial envelope to
                    // their requested weight before any admissibility decision.
                    auto sensed = corrected.Evaluate(sample,m_switch,AWFFailurePolicy::UseOriginal,p.wallRe >= m_minimumRe);
                    const double requested = m_mask[q]*sensed.requestedWeight;
                    p.selected = corrected.EvaluateFixedWeight(th,force,requested,m_failurePolicy,p.wallRe >= m_minimumRe);
                    p.selected.departure = sensed.departure; p.selected.sensors = sensed.sensors;
                }
                else p.selected = corrected.EvaluateFixedWeight(th,force,m_mask[q]*m_fixedWeight,m_failurePolicy,p.wallRe >= m_minimumRe);
            }
            const auto metric = IntegrateAWFMetric(wall,th,force,frame.curvature,p.selected.weight,numerator,cutoff);
            p.p = metric.meanProduction; p.d = metric.meanDissipation;
            eps[q] = p.selected.values.dissipationH; p.epsH = eps[q];
            if (!(p.selected.values.dissipationP > m_minimumEpsilon)) throw std::invalid_argument("terrain epsilon P guard");
            tangentMass[q] = alpha*V/(J*(1-w));
            for (int i = 0; i < 2; ++i)
            {
                const double pred = m_strip.EliminateVelocityP(bth[i],btn[i]);
                tangentRhs[q][i] = alpha*V/J*(w*tn[i]/(1-w)+pred)-V/J*pressure[i]-p.selected.values.traction[i]/J;
            }
            // Projection of the auxiliary normal-shear flux in the rotating
            // streamwise frame: integral c*T_tn d eta. Span has zero curvature.
            tangentRhs[q][0] -= frame.curvature*metric.integratedShear[0]/J;
            kMass[q] = alpha*V/(J*(1-wk));
            kRhs[q] = alpha*V/J*(wk*n[m_dim][q]/(1-wk)+m_strip.EliminateKP(bh[m_dim][q],bn[m_dim][q]))+V/J*(p.p-p.d);
            p.tangentSource = J*tangentRhs[q][0]; p.kSource = J*kRhs[q];
        }
    }
    catch (const std::exception &) { valid = false; }
    Require(valid,"nonpositive eliminated state or invalid AWF/metric reference; no floors or cutoff clipping");
    m_projectedError = ProjectConstraint(m_normalError);
    const auto projectedTarget = AdjointConstraint(ProjectConstraint(m_normalTarget));
    // Augment with C* C and load the reaction through C*. This preserves
    // symmetry/positivity and makes augmentation vanish at the same
    // conservative constraint used by the multiplier update.
    m_beta = 0;
    for (double value : beta) m_beta = std::max(m_beta,value);
    m_session->GetComm()->AllReduce(m_beta,LibUtilities::ReduceMax);
    double constraint = 0;
    for (int q = 0; q < m_localPoints; ++q)
        m_multiplier[q] -= m_normalRelaxation*m_beta*m_projectedError[q];
    const auto reaction = AdjointConstraint(m_multiplier);
    for (int q = 0; q < m_localPoints; ++q)
    {
        constraint = std::max(constraint,std::abs(m_projectedError[q]));
        const auto &normal = m_frames[q].normal;
        const auto g = Cartesian(tangentRhs[q],q);
        m_velocityMass[q] = tangentMass[q];
        for (int i = 0; i < m_dim; ++i)
        {
            const double ni = Component(normal,i);
            m_diagonal[i][q] = tangentMass[q];
            m_rhs[i][q] = Component(g,i)+(m_beta*projectedTarget[q]+reaction[q])*ni;
        }
        m_diagonal[m_dim][q] = kMass[q]; m_rhs[m_dim][q] = kRhs[q];
    }
    m_session->GetComm()->AllReduce(constraint,LibUtilities::ReduceMax); m_constraint = constraint;
    for (int i = 0; i <= m_dim; ++i)
    {
        if (!m_log || i < m_dim) SetBoundary(i,m_rhs[i],false);
        // A constant positive diagonal is used only by the cached native
        // preconditioner. The matrix-free operator applies the actual diagonal.
        for (int b : m_boundaries)
            if (m_fields[i]->GetBndConditions()[b]->GetBoundaryConditionType() == SpatialDomains::eRobin)
            {
                auto bc = std::static_pointer_cast<SpatialDomains::RobinBoundaryCondition>(m_fields[i]->GetBndConditions()[b]);
                bc->m_robinPrimitiveCoeff = std::make_shared<LibUtilities::Equation>(m_session->GetInterpreter(),"1");
            }
    }
    SetBoundary(m_dim+1,eps,true);
}

void RANSTerrainWall::AddRobinAuditResidual(int f, const Fields &state, Array<OneD,NekDouble> &residual) const
{
    if (f > m_dim || (m_log && f == m_dim)) return;
    if (f < m_dim)
    {
        if (f == 0 || m_robinNext != f)
        {
            Fields volume(m_dim);
            for (int i = 0; i < m_dim; ++i) volume[i] = state[i];
            m_robinVelocity = Trace(volume);
            m_robinNormal = Array<OneD,NekDouble>(m_localPoints,0.0);
            for (int q = 0; q < m_localPoints; ++q)
                for (int i = 0; i < m_dim; ++i)
                    m_robinNormal[q] += Component(m_frames[q].normal,i)*m_robinVelocity[i][q];
            m_robinProjected = AdjointConstraint(ProjectConstraint(m_robinNormal));
        }
        m_robinNext = f+1;
    }
    const auto map = std::dynamic_pointer_cast<MultiRegions::ContField>(m_fields[f])->GetLocalToGlobalMap();
    int offset = 0, wall = 0;
    for (int b = 0; b < m_fields[f]->GetBndCondExpansions().size(); ++b)
    {
        auto exp = m_fields[f]->GetBndCondExpansions()[b];
        if (m_fields[f]->GetBndConditions()[b]->GetUserDefined() == m_tag)
        {
            Array<OneD,NekDouble> trace, coefficients(exp->GetNcoeffs());
            m_fields[f]->ExtractPhysToBnd(b,state[f],trace);
            const int start = m_offsets[wall++];
            if (f < m_dim)
                for (int q = 0; q < trace.size(); ++q)
                {
                    const int k = start+q;
                    const double ni = Component(m_frames[k].normal,f);
                    trace[q] = m_velocityMass[k]*(trace[q]-ni*m_robinNormal[k])+m_beta*ni*m_robinProjected[k];
                }
            else Vmath::Vmul(trace.size(),trace,1,m_diagonal[f]+start,1,trace,1);
            exp->IProductWRTBase(trace,coefficients);
            for (int j = 0; j < coefficients.size(); ++j)
                residual[map->GetBndCondCoeffsToLocalCoeffsMap()[offset+j]] += map->GetBndCondCoeffsToLocalCoeffsSign()[offset+j]*coefficients[j];
        }
        offset += exp->GetNcoeffs();
    }
}

void RANSTerrainWall::WritePointDiagnostics(NekDouble time, int phase)
{
    if (!m_pointDiagnostics.is_open()) return;
    for (int q = 0; q < m_localPoints; ++q)
    {
        const auto &f = m_frames[q]; const auto &p = m_points[q]; const auto &s = p.sample; const auto &v = p.selected.values;
        m_pointDiagnostics << time << ',' << phase << ',' << m_globalOffset+q << ',' << m_x[q] << ',' << m_span[q]
            << ',' << f.curvature << ',' << f.Jacobian(m_strip.height) << ',' << f.Volume(m_strip.height)
            << ',' << m_strip.height << ',' << m_strip.pointP << ',' << m_strip.pointN << ',' << m_roughness[q]
            << ',' << m_logRoughness[q] << ',' << m_mask[q] << ',' << p.kP << ',' << s.matchingVelocity[0]
            << ',' << s.velocityP[0] << ',' << s.velocityN[0] << ',' << s.kN << ',' << s.epsilonN
            << ',' << s.pressureForcing[0]+s.convectionForcing[0] << ',' << s.pressureForcing[1]+s.convectionForcing[1]
            << ',' << v.traction[0] << ',' << v.traction[1] << ',' << p.p << ',' << p.d
            << ',' << v.dissipationP << ',' << p.epsH << ',' << p.selected.weight << ',' << p.selected.requestedWeight
            << ',' << p.selected.departure << ',' << static_cast<int>(p.selected.reason) << ',' << p.wallRe << ',' << p.cutoff
            << ',' << m_normalTarget[q] << ',' << m_normalError[q] << ',' << m_projectedError[q] << ',' << m_multiplier[q]
            << ',' << p.tangentSource << ',' << p.kSource << ',' << s.matchingVelocity[1] << ',' << s.velocityP[1]
            << ',' << s.velocityN[1] << ',' << s.pressureForcing[0] << ',' << s.pressureForcing[1]
            << ',' << s.convectionForcing[0] << ',' << s.convectionForcing[1] << ',' << s.normalVelocityP;
        for (double value : p.selected.sensors) m_pointDiagnostics << ',' << value;
        m_pointDiagnostics << ',' << (m_equilibrium && m_switched && m_mask[q] > 0)
            << ',' << p.resolvedShear << ',' << p.tangentFluxResidual << ',' << p.kFluxResidual << '\n';
    }
    m_pointDiagnostics << std::flush;
}
void RANSTerrainWall::WriteDiagnostics(NekDouble time, const Fields &state, const Fields &, const Fields &,
    const Array<OneD,NekDouble> &)
{
    const auto h = Trace(state);
    Fields derivativeTraces((m_dim+1)*m_dim);
    for (int i = 0; i <= m_dim; ++i)
        for (int j = 0; j < m_dim; ++j)
        {
            Fields derivative(1); derivative[0] = Array<OneD,NekDouble>(m_fields[0]->GetTotPoints());
            m_fields[i]->PhysDeriv(DirCartesianMap[j],state[i],derivative[0]);
            derivativeTraces[i*m_dim+j] = Trace(derivative)[0];
        }
    std::string diffusion; m_session->LoadSolverInfo("RANSResolvedDiffusion",diffusion,"Molecular");
    const double nu = diffusion == "Molecular" ? m_awf.molecularNu : 0;
    bool valid = true;
    try
    {
        for (int q = 0; q < m_localPoints; ++q)
        {
            RANSTensor gradient{};
            for (int i = 0; i < m_dim; ++i)
                for (int j = 0; j < m_dim; ++j)
                    gradient[m_dim == 2 && i == 1 ? 2 : i][m_dim == 2 && j == 1 ? 2 : j] = derivativeTraces[i*m_dim+j][q];
            const auto closure = EvaluateRANSClosure(m_reference.Model(),h[m_dim][q],h[m_dim+1][q],
                gradient,m_strip.height,m_reference.Parameters(),nu);
            const auto &t = m_frames[q].tangent; const auto &n = m_frames[q].normal;
            double shear = 0, kFlux = 0;
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    shear += t[i]*(closure.momentumDiffusion*(gradient[i][j]+gradient[j][i])-
                        closure.quadraticStress[i][j]-closure.cubicStress[i][j])*n[j];
            for (int i = 0; i < m_dim; ++i)
                kFlux += closure.kDiffusion*Component(n,i)*derivativeTraces[m_dim*m_dim+i][q];
            RANSVector right{};
            for (int i = 0; i < m_dim; ++i) right[m_dim == 2 && i == 1 ? 2 : i] = m_rhs[i][q];
            m_points[q].resolvedShear = shear;
            m_points[q].tangentFluxResidual = -shear+m_velocityMass[q]*Dot(Vector(h,q),t)-Dot(right,t);
            m_points[q].kFluxResidual = -kFlux+m_diagonal[m_dim][q]*h[m_dim][q]-m_rhs[m_dim][q];
        }
    }
    catch (const std::exception &) { valid = false; }
    Require(valid,"cannot audit the resolved interface stress/energy flux");
    double pointwise = 0, regularPointwise = 0, projection = 0, chiMin = 1, chiMax = 0;
    double kMin = std::numeric_limits<double>::max(), eMin = kMin;
    Array<OneD,NekDouble> ones(m_localPoints,1), wallArea(m_localPoints);
    for (int q = 0; q < m_localPoints; ++q)
    {
        pointwise = std::max(pointwise,std::abs(m_normalError[q]));
        if (m_terrain.terrain != RANSTerrain::CosineHill ||
            std::abs(std::abs(m_x[q]-m_terrain.centre)-m_terrain.hillHalfLength) >= .05*m_terrain.hillHalfLength)
            regularPointwise = std::max(regularPointwise,std::abs(m_normalError[q]));
        projection = std::max(projection,std::abs(h[m_dim+1][q]-m_points[q].epsH));
        chiMin = std::min(chiMin,m_points[q].selected.weight); chiMax = std::max(chiMax,m_points[q].selected.weight);
        kMin = std::min(kMin,m_points[q].kP); eMin = std::min(eMin,h[m_dim+1][q]);
        wallArea[q] = 1/m_frames[q].Jacobian(m_strip.height);
    }
    for (auto p : {&pointwise,&regularPointwise,&projection,&chiMax}) m_session->GetComm()->AllReduce(*p,LibUtilities::ReduceMax);
    for (auto p : {&chiMin,&kMin,&eMin}) m_session->GetComm()->AllReduce(*p,LibUtilities::ReduceMin);
    const double area = SurfaceIntegral(ones), physical = SurfaceIntegral(wallArea);
    const double normalDefect = SurfaceIntegral(m_normalError);
    Array<OneD,NekDouble> div(m_fields[0]->GetTotPoints(),0.0), derivative(div.size());
    for (int i = 0; i < m_dim; ++i)
    {
        m_fields[i]->PhysDeriv(DirCartesianMap[i],state[i],derivative);
        Vmath::Vadd(div.size(),div,1,derivative,1,div,1);
    }
    const double mass = m_fields[0]->Integral(div);
    for (int q = 0; q < m_localPoints; ++q) wallArea[q] = m_points[q].tangentFluxResidual*m_points[q].tangentFluxResidual;
    const double tangentResidual = std::sqrt(std::max(0.0,SurfaceIntegral(wallArea))/area);
    for (int q = 0; q < m_localPoints; ++q) wallArea[q] = m_points[q].kFluxResidual*m_points[q].kFluxResidual;
    const double kResidual = std::sqrt(std::max(0.0,SurfaceIntegral(wallArea))/area);
    for (int q = 0; q < m_localPoints; ++q) wallArea[q] = m_normalError[q]*m_normalError[q];
    const double normalL2 = std::sqrt(std::max(0.0,SurfaceIntegral(wallArea))/area);
    if (m_diagnostics.is_open())
        m_diagnostics << time << ',' << m_constraint << ',' << pointwise << ',' << projection << ',' << m_globalPoints
            << ',' << m_remoteSamples << ',' << chiMin << ',' << chiMax << ',' << kMin << ',' << eMin
            << ',' << area << ',' << physical << ',' << m_samplingResidual << ',' << m_projectionAudit
            << ',' << m_resolvedVolume << ',' << m_unresolvedVolume << ',' << normalDefect << ',' << mass
            << ',' << tangentResidual << ',' << kResidual << ',' << normalL2 << ',' << regularPointwise << '\n' << std::flush;
}
}
