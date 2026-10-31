///////////////////////////////////////////////////////////////////////////////
// File: RANSAnalyticalWallStrip.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSAnalyticalWallStrip.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/ContField.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace Nektar
{
std::string RANSAnalyticalWallStrip::RestartSignature() const
{
    std::ostringstream s;
    s << std::setprecision(17) << m_wallModel << ':' << m_geometry.height << ':'
      << m_geometry.pointP << ':' << m_geometry.pointN << ':' << m_awfGeometry.molecularNu
      << ':' << m_awfGeometry.roughness << ':' << m_parameters.cMu << ':'
      << m_parameters.cLength << ':' << m_parameters.smoothIntercept << ':'
      << m_minimumK << ':' << m_minimumEpsilon;
    if (m_equilibrium)
    {
        s << ':' << m_switched << ':' << m_fixedWeight << ':' << m_minimumRe
          << ':' << static_cast<int>(m_failurePolicy);
        if (m_switched)
        {
            s << ':' << m_switch.on << ':' << m_switch.off;
            for (auto v : m_switch.weights) s << ':' << v;
            for (auto v : m_switch.incomingDirection) s << ':' << v;
        }
    }
    if (m_spatial.enabled) s << ":spatial:" << m_spatial.begin << ':' << m_spatial.end;
    s << RANSWallPropertySignature(m_session);
    return s.str();
}
using namespace SolverUtils;
using MultiRegions::DirCartesianMap;

void RANSAnalyticalWallStrip::Require(bool valid, const std::string &message) const
{
    int bad = valid ? 0 : 1;
    m_session->GetComm()->AllReduce(bad, LibUtilities::ReduceMax);
    ASSERTL0(!bad, "RANS AWF strip: " + message);
}

RANSAnalyticalWallStrip::RANSAnalyticalWallStrip(
    const LibUtilities::SessionReaderSharedPtr &session,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    int dimension, NekDouble wallOrigin, const RANSNeutralReference &reference)
    : m_session(session), m_fields(fields), m_reference(reference),
      m_dim(dimension), m_localPoints(0),
      m_globalPoints(0), m_globalOffset(0), m_remoteSamples(0)
{
    m_session->LoadSolverInfo("RANSWallModel", m_wallModel, "AWFOriginal");
    m_equilibrium = m_wallModel == "AWFEquilibrium";
    m_spatial = ReadRANSWallTransition(m_session,ReadRANSWallGeometry(m_session,wallOrigin));
    ASSERTL0(!m_spatial.enabled && !m_session->DefinesFunction("RANSWallProperties","roughness") &&
                 !m_session->DefinesFunction("RANSWallProperties","z0"),
             "Spatial AWF transitions/properties require RANSWallGeometry=Plane or CosineHill");
    RANSWallPropertySignature(m_session); // validate stationary expressions
    m_failurePolicy = AWFFailurePolicy::UseOriginal;
    if (m_equilibrium)
    {
        ASSERTL0(!m_session->DefinesParameter("RANSZ0") && !m_session->DefinesFunction("RANSWallProperties","z0"),
            "AWFEquilibrium derives aerodynamic roughness from its actual matching anchor; do not prescribe RANSZ0");
        std::string weightMode, failure;
        m_session->LoadSolverInfo("RANSAWFWeight", weightMode, "Switched");
        m_session->LoadSolverInfo("RANSAWFFailurePolicy", failure, "UseOriginal");
        ASSERTL0(weightMode == "Switched" || weightMode == "Fixed", "RANSAWFWeight must be Switched or Fixed");
        ASSERTL0(failure == "UseOriginal" || failure == "RequirePreservation", "Invalid RANSAWFFailurePolicy");
        m_switched = weightMode == "Switched";
        m_failurePolicy = failure == "RequirePreservation" ? AWFFailurePolicy::RequirePreservation : AWFFailurePolicy::UseOriginal;
        m_session->LoadParameter("RANSAWFFixedWeight", m_fixedWeight, 1.0);
        m_session->LoadParameter("RANSAWFMinimumRe", m_minimumRe, 100.0);
        ASSERTL0(std::isfinite(m_fixedWeight) && m_fixedWeight >= 0 && m_fixedWeight <= 1 &&
                 std::isfinite(m_minimumRe) && m_minimumRe > 0, "Invalid AWF weight or high-Re threshold");
        if (m_switched)
        {
            // The note does not calibrate thresholds. Require them explicitly.
            m_session->LoadParameter("RANSAWFOn", m_switch.on);
            m_session->LoadParameter("RANSAWFOff", m_switch.off);
            const std::array<std::string, 7> names{{"U", "DeltaU", "K", "Epsilon", "Pressure", "Convection", "Normal"}};
            for (int i = 0; i < names.size(); ++i)
                m_session->LoadParameter("RANSAWFWeight" + names[i], m_switch.weights[i], 1.0);
            m_session->LoadParameter("RANSAWFDirectionX", m_switch.incomingDirection[0], 1.0);
            m_session->LoadParameter("RANSAWFDirectionY", m_switch.incomingDirection[1], 0.0);
            m_switch.Validate();
            ASSERTL0(m_dim != 2 || m_switch.incomingDirection[1] == 0,
                "A 2D x-z wall requires a streamwise incoming direction");
        }
    }
    int pointDiagnostics = 0;
    m_session->LoadParameter("RANSAWFPointDiagnostics", pointDiagnostics, 0);
    ASSERTL0(pointDiagnostics == 0 || pointDiagnostics == 1, "RANSAWFPointDiagnostics must be 0 or 1");
    m_writePoints = pointDiagnostics;
    m_session->LoadParameter("RANSStripHeight", m_geometry.height, 2.0);
    m_session->LoadParameter("RANSStripPointP", m_geometry.pointP, 1.0);
    m_session->LoadParameter("RANSStripPointN", m_geometry.pointN, 4.0);
    m_session->LoadParameter("Kinvis", m_awfGeometry.molecularNu);
    m_session->LoadParameter("RANSRoughnessHeight", m_awfGeometry.roughness, .1);
    m_session->LoadParameter("RANSAWFCmu", m_parameters.cMu, .09);
    m_session->LoadParameter("RANSAWFCLength", m_parameters.cLength, 2.55);
    m_session->LoadParameter("RANSAWFSmoothIntercept", m_parameters.smoothIntercept, 10.7);
    m_session->LoadParameter("RANSMinimumK", m_minimumK, 0.0);
    m_session->LoadParameter("RANSMinimumEpsilon", m_minimumEpsilon, 0.0);
    m_awfGeometry.height = m_geometry.height;
    m_awfGeometry.pointP = m_geometry.pointP;
    m_geometry.Validate(); m_awfGeometry.Validate(); m_parameters.Validate();
    bool valid = true;
    for (int e = 0; e < m_fields[0]->GetExpSize(); ++e)
        valid &= m_fields[0]->GetExp(e)->DetShapeType() ==
            (m_dim == 2 ? LibUtilities::eQuadrilateral : LibUtilities::eHexahedron);
    Require(valid, "the initial flat-strip adapter requires quadrilaterals or hexahedra");
    for (int b = 0; b < m_fields[0]->GetBndConditions().size(); ++b)
        if (m_fields[0]->GetBndConditions()[b]->GetUserDefined() == m_wallModel)
        {
            m_boundaries.push_back(b); m_offsets.push_back(m_localPoints);
            auto boundary = m_fields[0]->GetBndCondExpansions()[b];
            const int count = boundary->GetTotPoints();
            Array<OneD, NekDouble> x(count, 0.0), y(count, 0.0), z(count, 0.0);
            Fields normals;
            boundary->GetCoords(x, y, z);
            m_fields[0]->GetBoundaryNormals(b, normals);
            for (int q = 0; q < count; ++q)
                valid &= std::abs((m_dim == 2 ? y[q] : z[q]) - wallOrigin - m_geometry.height) < 1e-9 &&
                         std::abs(normals[m_dim - 1][q] + 1) < 1e-10;
            m_localPoints += count;
        }
    Require(valid, "the analytical wall must occupy the flat lower interface at wallOrigin+H");
    m_wallX = Array<OneD,NekDouble>(m_localPoints);
    m_roughness = Array<OneD,NekDouble>(m_localPoints,m_awfGeometry.roughness);
    m_mask = Array<OneD,NekDouble>(m_localPoints,1.0);
    for (int j = 0; j < m_boundaries.size(); ++j)
    {
        auto boundary = fields[0]->GetBndCondExpansions()[m_boundaries[j]];
        const int count = boundary->GetTotPoints(), offset = m_offsets[j];
        Array<OneD,NekDouble> x(count,0.0),y(count,0.0),z(count,0.0),property(count);
        boundary->GetCoords(x,y,z);
        for (int q = 0; q < count; ++q)
        {
            m_wallX[offset+q] = x[q]; m_mask[offset+q] = m_spatial.Weight(x[q]);
            if (m_dim == 2) y[q] = wallOrigin; else z[q] = wallOrigin;
        }
        if (m_session->DefinesFunction("RANSWallProperties","roughness"))
        {
            m_session->GetFunction("RANSWallProperties","roughness")->Evaluate(x,y,z,0,property);
            Vmath::Vcopy(count,property.data(),1,m_roughness.data()+offset,1);
        }
    }
    try
    {
        for (double roughness : m_roughness)
        {
            auto g = m_awfGeometry; g.roughness = roughness; g.Validate();
        }
    }
    catch (const std::exception &) { valid = false; }
    Require(valid,"invalid spatial AWF roughness");
    auto comm = m_session->GetComm();
    Array<OneD, int> localCount(1, m_localPoints), counts, offsets(comm->GetSize());
    comm->AllGather(localCount, counts);
    for (int r = 0; r < counts.size(); ++r)
    {
        offsets[r] = m_globalPoints;
        m_wallRanks += counts[r] > 0;
        if (r == comm->GetRank()) m_globalOffset = m_globalPoints;
        m_globalPoints += counts[r];
    }
    Require(m_globalPoints > 0, "the analytical wall requires a matching tagged wall region");
    Fields coordinates(3), searchShift(3);
    for (int d = 0; d < 3; ++d)
    {
        Array<OneD, NekDouble> local(m_localPoints, 0.0), shift(m_localPoints, 0.0);
        for (int j = 0; j < m_boundaries.size(); ++j)
        {
            auto boundary = m_fields[0]->GetBndCondExpansions()[m_boundaries[j]];
            const int count = boundary->GetTotPoints();
            Array<OneD, NekDouble> x(count, 0.0), y(count, 0.0), z(count, 0.0);
            boundary->GetCoords(x, y, z);
            const auto &data = d == 0 ? x : (d == 1 ? y : z);
            Vmath::Vcopy(count, data.data(), 1, local.data() + m_offsets[j], 1);
            if (d < m_dim - 1)
                for (int e = 0; e < boundary->GetExpSize(); ++e)
                {
                    const int offset = boundary->GetPhys_Offset(e);
                    const int size = boundary->GetExp(e)->GetTotPoints();
                    const double centre = .5 *
                        (*std::min_element(data.begin() + offset, data.begin() + offset + size) +
                         *std::max_element(data.begin() + offset, data.begin() + offset + size));
                    for (int q = offset; q < offset + size; ++q)
                        shift[m_offsets[j] + q] = 1e-8 * (centre - data[q]);
                }
        }
        coordinates[d] = Array<OneD, NekDouble>(m_globalPoints);
        comm->AllGatherv(local, coordinates[d], counts, offsets);
        searchShift[d] = Array<OneD, NekDouble>(m_globalPoints);
        comm->AllGatherv(shift, searchShift[d], counts, offsets);
    }
    Array<OneD, int> ids(m_globalPoints, std::numeric_limits<int>::max());
    m_sampleElements.resize(m_globalPoints, -1);
    m_sampleCoordinates.resize(m_globalPoints);
    // Select the element just ABOVE N, then evaluate at the exact N height.
    // A search shifted into the wall patch selects its matching column at
    // shared tangential vertices. This preserves one-sided derivatives there;
    // choosing the same neighbour for BOTH patches would spoil strip continuity.
    // Only ownership uses the shift: evaluation remains at the exact N point.
    for (int q = 0; q < m_globalPoints; ++q)
    {
        Array<OneD, NekDouble> point(m_dim), local(m_dim);
        for (int d = 0; d < m_dim; ++d) point[d] = coordinates[d][q] + searchShift[d][q];
        point[m_dim - 1] = wallOrigin + m_geometry.pointN + 1e-8 * m_geometry.height;
        const int e = m_fields[0]->GetExpIndex(point, local, 1e-11);
        if (e >= 0)
        {
            for (int d = 0; d < m_dim; ++d) point[d] = coordinates[d][q];
            point[m_dim - 1] = wallOrigin + m_geometry.pointN;
            m_fields[0]->GetExp(e)->GetGeom()->GetLocCoords(point, local);
            m_sampleElements[q] = e;
            m_sampleCoordinates[q] = local;
            ids[q] = m_fields[0]->GetExp(e)->GetGeom()->GetGlobalID();
        }
    }
    comm->AllReduce(ids, LibUtilities::ReduceMin);
    Array<OneD, int> owners(m_globalPoints, std::numeric_limits<int>::max());
    valid = true;
    for (int q = 0; q < m_globalPoints; ++q)
    {
        valid &= ids[q] != std::numeric_limits<int>::max();
        const int e = m_sampleElements[q];
        if (e >= 0 && m_fields[0]->GetExp(e)->GetGeom()->GetGlobalID() == ids[q])
            owners[q] = comm->GetRank();
        else m_sampleElements[q] = -1;
    }
    Require(valid, "every physical N point must be inside the CG domain, below its top");
    comm->AllReduce(owners, LibUtilities::ReduceMin);
    for (int q = 0; q < m_localPoints; ++q)
    {
        const int owner = owners[m_globalOffset + q];
        m_remoteSamples += owner != comm->GetRank();
        m_remoteWithoutWallSamples += counts[owner] == 0;
    }
    comm->AllReduce(m_remoteSamples, LibUtilities::ReduceSum);
    comm->AllReduce(m_remoteWithoutWallSamples, LibUtilities::ReduceSum);
    for (auto p : {&m_production, &m_dissipation, &m_epsilonP, &m_epsilonH,
                   &m_kP, &m_tractionX, &m_tractionY, &m_forcingX})
        *p = Array<OneD, NekDouble>(m_localPoints, 0.0);
    m_points.resize(m_localPoints);
    if (m_writePoints)
    {
        m_pointDiagnostics.open(m_session->GetSessionName() + ".awf-points-r" + std::to_string(comm->GetRank()) + ".csv");
        ASSERTL0(m_pointDiagnostics.good(), "Cannot open AWF point diagnostics");
        m_pointDiagnostics << std::setprecision(17)
            << "time,phase,point,kP,uH,vH,uP,vP,uN,vN,kN,epsN,wP,pressure_x,pressure_y,convection_x,convection_y,forcing_x,forcing_y,traction_x,traction_y,production,dissipation,epsP,epsH,original_dissipation,original_epsP,original_epsH,reference_u,reference_A,reference_P0,lambda,minimum_lambda,cutoff,log_z0,wall_Re,admissible,requested_chi,chi,departure,reason,r_U,r_DeltaU,r_K,r_Epsilon,r_Pressure,r_Convection,r_Normal,x,mask,roughness\n";
    }
    if (comm->GetRank() == 0)
    {
        m_diagnostics.open(m_session->GetSessionName() + ".awf.csv");
        ASSERTL0(m_diagnostics.good(), "Cannot open AWF strip diagnostics");
        m_diagnostics << std::setprecision(17)
            << "time,area,traction_x,traction_y,mean_production,mean_dissipation,min_kP,min_epsP,epsilonH_error,mass_error,normal_velocity_max,momentum_budget_x,momentum_budget_y,k_budget,budget_available,global_wall_points,remote_N_samples,chi,guard_activations,mean_uH,mean_uN,mean_uP,mean_kP,mean_forcing_x,mean_epsP,mean_epsH,wall_ranks,remote_N_without_wall_samples,epsilonH_projection_error,chi_min,chi_max,requested_chi,departure_max,lambda_min,cutoff_max,wall_Re_min,reference_P0,original_dissipation,active_points,inactive_points,cutoff_fallback_points,regime_fallback_points,range_fallback_points\n";
    }
}

RANSAnalyticalWallStrip::Fields RANSAnalyticalWallStrip::Trace(const Fields &values) const
{
    Fields result(values.size());
    for (int f = 0; f < values.size(); ++f)
    {
        result[f] = Array<OneD, NekDouble>(m_localPoints);
        for (int j = 0; j < m_boundaries.size(); ++j)
        {
            Array<OneD, NekDouble> trace;
            m_fields[f < m_fields.size() ? f : 0]->ExtractPhysToBnd(m_boundaries[j], values[f], trace);
            Vmath::Vcopy(trace.size(), trace.data(), 1, result[f].data() + m_offsets[j], 1);
        }
    }
    return result;
}

RANSAnalyticalWallStrip::Fields RANSAnalyticalWallStrip::SampleN(const Fields &values) const
{
    Fields result(values.size());
    for (int f = 0; f < values.size(); ++f)
    {
        Array<OneD, NekDouble> all(m_globalPoints, 0.0);
        for (int q = 0; q < m_globalPoints; ++q)
        {
            const int e = m_sampleElements[q];
            if (e >= 0)
                all[q] = m_fields[0]->GetExp(e)->StdPhysEvaluate(
                    m_sampleCoordinates[q], values[f] + m_fields[0]->GetPhys_Offset(e));
        }
        m_session->GetComm()->AllReduce(all, LibUtilities::ReduceSum);
        result[f] = Array<OneD, NekDouble>(m_localPoints);
        Vmath::Vcopy(m_localPoints, all + m_globalOffset, 1, result[f], 1);
    }
    return result;
}

void RANSAnalyticalWallStrip::SetBoundary(int field,
    const Array<OneD, NekDouble> &values, bool dirichlet)
{
    for (int j = 0; j < m_boundaries.size(); ++j)
    {
        auto boundary = m_fields[field]->GetBndCondExpansions()[m_boundaries[j]];
        Array<OneD, NekDouble> local = values + m_offsets[j];
        if (dirichlet) boundary->FwdTransBndConstrained(local, boundary->UpdateCoeffs());
        else boundary->IProductWRTBase(local, boundary->UpdateCoeffs());
    }
}

void RANSAnalyticalWallStrip::CheckInitialState(const Fields &state) const
{
    if (!m_equilibrium || m_failurePolicy != AWFFailurePolicy::RequirePreservation ||
        (!m_switched && m_fixedWeight == 0)) return;
    const auto h = Trace(state), n = SampleN(state);
    bool valid = true;
    try
    {
        for (int q = 0; q < m_localPoints; ++q)
        {
            if (m_mask[q] == 0) continue;
            const double k = m_geometry.EliminateKP(h[m_dim][q], n[m_dim][q]);
            auto g = m_awfGeometry; g.roughness = m_roughness[q];
            EquilibriumWallFunction wall(g, k, m_reference, m_parameters);
            valid &= wall.Reference().admissible && std::sqrt(k) * m_geometry.height /
                m_awfGeometry.molecularNu >= m_minimumRe;
        }
    }
    catch (const std::exception &) { valid = false; }
    Require(valid, "requested equilibrium reference fails cutoff/high-Re preflight; no cutoff clipping is allowed");
}

void RANSAnalyticalWallStrip::AddRobinAuditResidual(int field,
    const Fields &state, Array<OneD, NekDouble> &residual) const
{
    if (field > m_dim || field == m_dim - 1) return;
    const double weight = field == m_dim ? m_geometry.KWeightN() : m_geometry.VelocityWeightN();
    const double coefficient = m_geometry.height * m_alpha / (1 - weight);
    const auto map = std::dynamic_pointer_cast<MultiRegions::ContField>(m_fields[field])->GetLocalToGlobalMap();
    int offset = 0;
    for (int b = 0; b < m_fields[field]->GetBndCondExpansions().size(); ++b)
    {
        auto boundary = m_fields[field]->GetBndCondExpansions()[b];
        if (m_fields[field]->GetBndConditions()[b]->GetUserDefined() == m_wallModel)
        {
            Array<OneD, NekDouble> trace, coefficients(boundary->GetNcoeffs());
            m_fields[field]->ExtractPhysToBnd(b, state[field], trace);
            Vmath::Smul(trace.size(), coefficient, trace, 1, trace, 1);
            boundary->IProductWRTBase(trace, coefficients);
            for (int j = 0; j < coefficients.size(); ++j)
                residual[map->GetBndCondCoeffsToLocalCoeffsMap()[offset+j]] +=
                    map->GetBndCondCoeffsToLocalCoeffsSign()[offset+j] * coefficients[j];
        }
        offset += boundary->GetNcoeffs();
    }
}

void RANSAnalyticalWallStrip::WritePointDiagnostics(NekDouble time, int phase)
{
    if (!m_pointDiagnostics.is_open()) return;
    for (int q = 0; q < m_localPoints; ++q)
    {
        const auto &p = m_points[q]; const auto &s = p.samples;
        const auto &v = p.selected.values; const auto &r = p.reference;
        m_pointDiagnostics << time << ',' << phase << ',' << m_globalOffset + q << ',' << m_kP[q];
        for (const auto &vector : {s.matchingVelocity, s.velocityP, s.velocityN})
            for (double x : vector) m_pointDiagnostics << ',' << x;
        m_pointDiagnostics << ',' << s.kN << ',' << s.epsilonN << ',' << s.normalVelocityP;
        for (const auto &vector : {s.pressureForcing, s.convectionForcing})
            for (double x : vector) m_pointDiagnostics << ',' << x;
        for (int i = 0; i < 2; ++i) m_pointDiagnostics << ',' << s.pressureForcing[i] + s.convectionForcing[i];
        for (double x : v.traction) m_pointDiagnostics << ',' << x;
        for (double x : {v.meanProduction, v.meanDissipation, v.dissipationP, v.dissipationH,
                        p.original.meanDissipation, p.original.dissipationP, p.original.dissipationH,
                        r.matchingVelocity, r.numerator, r.production, r.lambda, r.minimumLambda,
                        r.cutoff, r.logAerodynamicRoughness, p.wallRe}) m_pointDiagnostics << ',' << x;
        m_pointDiagnostics << ',' << r.admissible << ',' << p.selected.requestedWeight << ','
            << p.selected.weight << ',' << p.selected.departure << ',' << static_cast<int>(p.selected.reason);
        for (double x : p.selected.sensors) m_pointDiagnostics << ',' << x;
        m_pointDiagnostics << ',' << m_wallX[q] << ',' << m_mask[q] << ',' << m_roughness[q];
        m_pointDiagnostics << '\n';
    }
    m_pointDiagnostics << std::flush;
}

void RANSAnalyticalWallStrip::Update(const Fields &state, const Fields &predictor,
    const Gradients &gradient, const Fields &pressureForce,
    const Fields &convection, const Fields &bodyForce, NekDouble alpha)
{
    m_volumePredictor = predictor;
    const auto h = Trace(state), n = SampleN(state);
    const auto bh = Trace(predictor), bn = SampleN(predictor);
    const auto ph = Trace(pressureForce), pn = SampleN(pressureForce);
    const auto ah = Trace(convection), an = SampleN(convection);
    const auto fh = Trace(bodyForce), fn = SampleN(bodyForce);
    m_predictorP = Fields(m_dim + 1); m_g = Fields(m_dim + 1); m_pressureP = Fields(m_dim - 1);
    for (int i = 0; i <= m_dim; ++i)
    {
        m_predictorP[i] = Array<OneD, NekDouble>(m_localPoints);
        m_g[i] = Array<OneD, NekDouble>(m_localPoints);
        if (i < m_dim - 1) m_pressureP[i] = Array<OneD, NekDouble>(m_localPoints);
    }
    Array<OneD, NekDouble> normal(m_localPoints, 0.0);
    for (int i = 0; i < m_dim - 1; ++i)
    {
        Fields derivative(1); derivative[0] = gradient[i][i];
        const auto dh = Trace(derivative), dn = SampleN(derivative);
        for (int q = 0; q < m_localPoints; ++q)
            normal[q] -= m_geometry.height * m_geometry.EliminateVelocityP(dh[0][q], dn[0][q]);
    }
    bool valid = true;
    try
    {
        for (int q = 0; q < m_localPoints; ++q)
        {
            m_kP[q] = m_geometry.EliminateKP(h[m_dim][q], n[m_dim][q]);
            if (!(m_kP[q] > m_minimumK)) throw std::invalid_argument("strip KP guard");
            auto geometry = m_awfGeometry; geometry.roughness = m_roughness[q];
            AnalyticalWallFunction wall(geometry, m_kP[q], m_parameters);
            auto &point = m_points[q];
            point = PointData{};
            auto &sample = point.samples;
            sample.heightN = m_geometry.pointN;
            sample.kN = n[m_dim][q]; sample.epsilonN = n[m_dim + 1][q];
            // The strip's constant tangential-divergence approximation gives
            // W(z)=z/H W_H, with impermeability at the physical wall.
            sample.normalVelocityP = m_geometry.pointP / m_geometry.height * h[m_dim - 1][q];
            AWFVector velocity{}, forcing{};
            for (int i = 0; i < m_dim - 1; ++i)
            {
                velocity[i] = h[i][q];
                m_pressureP[i][q] = m_geometry.EliminateVelocityP(ph[i][q], pn[i][q]);
                sample.matchingVelocity[i] = h[i][q];
                sample.velocityP[i] = m_geometry.EliminateVelocityP(h[i][q], n[i][q]);
                sample.velocityN[i] = n[i][q];
                sample.pressureForcing[i] = m_pressureP[i][q] - m_geometry.EliminateVelocityP(fh[i][q], fn[i][q]);
                sample.convectionForcing[i] = m_geometry.EliminateVelocityP(ah[i][q], an[i][q]);
                forcing[i] = sample.pressureForcing[i] + sample.convectionForcing[i];
            }
            point.original = wall.Evaluate(velocity, forcing);
            point.selected.values = point.original;
            point.selected.reason = AWFFallbackReason::SwitchInactive;
            if (m_equilibrium && m_mask[q] > 0)
            {
                EquilibriumWallFunction corrected(geometry, m_kP[q], m_reference, m_parameters);
                point.reference = corrected.Reference();
                point.wallRe = std::sqrt(m_kP[q]) * m_geometry.height / m_awfGeometry.molecularNu;
                const bool highRe = point.wallRe >= m_minimumRe;
                if (m_switched)
                {
                    const auto sensed = corrected.Evaluate(sample,m_switch,AWFFailurePolicy::UseOriginal,highRe);
                    point.selected = corrected.EvaluateFixedWeight(velocity,forcing,m_mask[q]*sensed.requestedWeight,m_failurePolicy,highRe);
                    point.selected.departure = sensed.departure; point.selected.sensors = sensed.sensors;
                }
                else point.selected = corrected.EvaluateFixedWeight(velocity,forcing,m_mask[q]*m_fixedWeight,m_failurePolicy,highRe);
            }
            const auto &values = point.selected.values;
            if (!(values.dissipationP > m_minimumEpsilon)) throw std::invalid_argument("strip EP guard");
            m_production[q] = values.meanProduction; m_dissipation[q] = values.meanDissipation;
            m_epsilonP[q] = values.dissipationP; m_epsilonH[q] = values.dissipationH;
            m_tractionX[q] = values.traction[0]; m_tractionY[q] = values.traction[1];
            m_forcingX[q] = forcing[0];
            for (int i = 0; i <= m_dim; ++i)
            {
                if (i == m_dim - 1) continue;
                const double weight = i == m_dim ? m_geometry.KWeightN() : m_geometry.VelocityWeightN();
                m_predictorP[i][q] = (bh[i][q] - weight * bn[i][q]) / (1 - weight);
                m_g[i][q] = m_geometry.height * alpha *
                    (weight * n[i][q] / (1 - weight) + m_predictorP[i][q]);
                if (i == m_dim)
                    m_g[i][q] += m_geometry.height * (values.meanProduction - values.meanDissipation);
                else m_g[i][q] -= m_geometry.height * m_pressureP[i][q] + values.traction[i];
            }
        }
    }
    catch (const std::exception &) { valid = false; }
    Require(valid, "nonpositive eliminated KP or inadmissible/nonfinite AWF state or requested equilibrium reference; no floor or cutoff clipping is applied");
    if (alpha != m_alpha)
        for (int b : m_boundaries)
            for (int i = 0; i <= m_dim; ++i)
            {
                if (i == m_dim - 1) continue;
                const double weight = i == m_dim ? m_geometry.KWeightN() : m_geometry.VelocityWeightN();
                std::ostringstream coefficient;
                coefficient << std::setprecision(17) << m_geometry.height * alpha / (1 - weight);
                auto bc = std::static_pointer_cast<SpatialDomains::RobinBoundaryCondition>(m_fields[i]->GetBndConditions()[b]);
                bc->m_robinPrimitiveCoeff = std::make_shared<LibUtilities::Equation>(
                    m_session->GetInterpreter(), coefficient.str());
            }
    m_alpha = alpha;
    for (int i = 0; i <= m_dim; ++i)
        if (i != m_dim - 1) SetBoundary(i, m_g[i], false);
    SetBoundary(m_dim - 1, normal, true);
    SetBoundary(m_dim + 1, m_epsilonH, true);
}

NekDouble RANSAnalyticalWallStrip::SurfaceIntegral(const Array<OneD, NekDouble> &values) const
{
    // Element integrals are local; the final collective includes wall-less ranks.
    double result = 0;
    for (int j = 0; j < m_boundaries.size(); ++j)
    {
        auto boundary = m_fields[0]->GetBndCondExpansions()[m_boundaries[j]];
        for (int e = 0; e < boundary->GetExpSize(); ++e)
            result += boundary->GetExp(e)->Integral(values + m_offsets[j] + boundary->GetPhys_Offset(e));
    }
    m_session->GetComm()->AllReduce(result, LibUtilities::ReduceSum);
    return result;
}

NekDouble RANSAnalyticalWallStrip::NaturalFlux(int field, bool &available) const
{
    double flux = 0;
    for (int b = 0; b < m_fields[field]->GetBndCondExpansions().size(); ++b)
    {
        const auto bc = m_fields[field]->GetBndConditions()[b];
        if (bc->GetBoundaryConditionType() == SpatialDomains::ePeriodic ||
            bc->GetUserDefined() == m_wallModel) continue;
        if (bc->GetBoundaryConditionType() != SpatialDomains::eNeumann)
        { available = false; continue; }
        auto boundary = m_fields[field]->GetBndCondExpansions()[b];
        Array<OneD, NekDouble> ones(boundary->GetTotPoints(), 1.0), coefficients(boundary->GetNcoeffs());
        boundary->FwdTransLocalElmt(ones, coefficients);
        flux += Vmath::Dot(coefficients.size(), coefficients, 1, boundary->GetCoeffs(), 1);
    }
    m_session->GetComm()->AllReduce(flux, LibUtilities::ReduceSum);
    return flux;
}

void RANSAnalyticalWallStrip::WriteDiagnostics(NekDouble time, const Fields &state,
    const Fields &gain, const Fields &loss, const Array<OneD, NekDouble> &pressure)
{
    const auto h = Trace(state), n = SampleN(state);
    Array<OneD, NekDouble> ones(m_localPoints, 1.0), work(m_localPoints), mass(m_localPoints, 0.0);
    const double area = SurfaceIntegral(ones);
    double epsError = 0, projectionError = 0, minK = std::numeric_limits<double>::max(), minE = minK, normalMax = 0;
    Array<OneD, NekDouble> imposedEpsilon(m_localPoints);
    for (int j = 0; j < m_boundaries.size(); ++j)
    {
        auto boundary = m_fields[m_dim + 1]->GetBndCondExpansions()[m_boundaries[j]];
        Array<OneD, NekDouble> local = imposedEpsilon + m_offsets[j];
        boundary->BwdTrans(boundary->GetCoeffs(), local);
    }
    for (int i = 0; i < m_dim - 1; ++i)
    {
        Fields derivative(1); derivative[0] = Array<OneD, NekDouble>(m_fields[0]->GetTotPoints());
        m_fields[0]->PhysDeriv(DirCartesianMap[i], state[i], derivative[0]);
        const auto dh = Trace(derivative), dn = SampleN(derivative);
        for (int q = 0; q < m_localPoints; ++q)
            mass[q] += m_geometry.height * m_geometry.EliminateVelocityP(dh[0][q], dn[0][q]);
    }
    double massError = 0;
    for (int q = 0; q < m_localPoints; ++q)
    {
        minK = std::min(minK, m_geometry.EliminateKP(h[m_dim][q], n[m_dim][q]));
        minE = std::min(minE, m_epsilonP[q]);
        epsError = std::max(epsError, std::abs(h[m_dim + 1][q] - imposedEpsilon[q]));
        projectionError = std::max(projectionError, std::abs(imposedEpsilon[q] - m_epsilonH[q]));
        massError = std::max(massError, std::abs(mass[q] + h[m_dim - 1][q]));
        normalMax = std::max(normalMax, std::abs(h[m_dim - 1][q]));
    }
    for (auto p : {&epsError, &projectionError, &massError, &normalMax}) m_session->GetComm()->AllReduce(*p, LibUtilities::ReduceMax);
    for (auto p : {&minK, &minE}) m_session->GetComm()->AllReduce(*p, LibUtilities::ReduceMin);
    // Constant-test budgets independently add outer volume rows and the strip
    // inventory/source row. Only report a budget if all other faces are natural
    // or periodic for that variable; prescribed-trace reactions are unknown.
    std::array<double, 3> budget{};
    bool available = true;
    const int nq = m_fields[0]->GetTotPoints();
    for (int i = 0; i <= m_dim; ++i)
    {
        if (i == m_dim - 1) continue;
        const double weight = i == m_dim ? m_geometry.KWeightN() : m_geometry.VelocityWeightN();
        for (int q = 0; q < m_localPoints; ++q)
        {
            const double p = (h[i][q] - weight * n[i][q]) / (1 - weight);
            work[q] = m_geometry.height * m_alpha * (p - m_predictorP[i][q]);
            work[q] += i == m_dim ? -m_geometry.height * (m_production[q] - m_dissipation[q]) :
                m_geometry.height * m_pressureP[i][q] + (i == 0 ? m_tractionX[q] : m_tractionY[q]);
        }
        double residual = SurfaceIntegral(work) - NaturalFlux(i, available);
        // The outer predictor contribution can be reconstructed from bH/bN
        // only on the boundary. Its full volume is retained below in Update.
        Array<OneD, NekDouble> volume(nq);
        for (int q = 0; q < nq; ++q) volume[q] = -gain[i][q] + loss[i][q] * state[i][q];
        if (i < m_dim)
        {
            Array<OneD, NekDouble> grad(nq);
            m_fields[m_fields.size() - 1]->PhysDeriv(DirCartesianMap[i], pressure, grad);
            Vmath::Vadd(nq, volume, 1, grad, 1, volume, 1);
        }
        // Add the full incoming-predictor inventory saved separately in Update.
        for (int q = 0; q < nq; ++q)
            volume[q] += m_alpha * (state[i][q] - m_volumePredictor[i][q]);
        residual += m_fields[0]->Integral(volume);
        budget[i == m_dim ? 2 : i] = residual;
    }
    int allAvailable = available ? 1 : 0;
    m_session->GetComm()->AllReduce(allAvailable, LibUtilities::ReduceMin);
    const double production = SurfaceIntegral(m_production) / area;
    const double dissipation = SurfaceIntegral(m_dissipation) / area;
    const double tx = SurfaceIntegral(m_tractionX), ty = SurfaceIntegral(m_tractionY);
    const double uh = SurfaceIntegral(h[0]) / area, un = SurfaceIntegral(n[0]) / area;
    for (int q = 0; q < m_localPoints; ++q)
        work[q] = m_geometry.EliminateVelocityP(h[0][q], n[0][q]);
    const double up = SurfaceIntegral(work) / area;
    const double kp = SurfaceIntegral(m_kP) / area, forcing = SurfaceIntegral(m_forcingX) / area;
    const double ep = SurfaceIntegral(m_epsilonP) / area, eh = SurfaceIntegral(m_epsilonH) / area;
    double chiMin = 1, chiMax = 0, departureMax = 0, cutoffMax = 0;
    double lambdaMin = std::numeric_limits<double>::max(), reMin = lambdaMin;
    Array<OneD, int> reasons(5, 0);
    Fields averaged(4);
    for (auto &v : averaged) v = Array<OneD, NekDouble>(m_localPoints);
    for (int q = 0; q < m_localPoints; ++q)
    {
        const auto &p = m_points[q];
        chiMin = std::min(chiMin, p.selected.weight); chiMax = std::max(chiMax, p.selected.weight);
        departureMax = std::max(departureMax, p.selected.departure);
        lambdaMin = std::min(lambdaMin, p.reference.lambda);
        cutoffMax = std::max(cutoffMax, p.reference.cutoff); reMin = std::min(reMin, p.wallRe);
        reasons[static_cast<int>(p.selected.reason)]++;
        averaged[0][q] = p.selected.weight; averaged[1][q] = p.selected.requestedWeight;
        averaged[2][q] = p.reference.production; averaged[3][q] = p.original.meanDissipation;
    }
    for (auto p : {&chiMin, &lambdaMin, &reMin}) m_session->GetComm()->AllReduce(*p, LibUtilities::ReduceMin);
    for (auto p : {&chiMax, &departureMax, &cutoffMax}) m_session->GetComm()->AllReduce(*p, LibUtilities::ReduceMax);
    m_session->GetComm()->AllReduce(reasons, LibUtilities::ReduceSum);
    std::array<double, 4> means;
    for (int i = 0; i < means.size(); ++i) means[i] = SurfaceIntegral(averaged[i]) / area;
    if (m_diagnostics.is_open())
        m_diagnostics << time << ',' << area << ',' << tx << ',' << ty << ',' << production << ',' << dissipation
            << ',' << minK << ',' << minE << ',' << epsError << ',' << massError << ',' << normalMax
            << ',' << budget[0] << ',' << budget[1] << ',' << budget[2] << ',' << allAvailable
            << ',' << m_globalPoints << ',' << m_remoteSamples << ',' << means[0] << ",0"
            << ',' << uh << ',' << un << ',' << up << ',' << kp << ',' << forcing << ',' << ep << ',' << eh
            << ',' << m_wallRanks << ',' << m_remoteWithoutWallSamples
            << ',' << projectionError
            << ',' << chiMin << ',' << chiMax << ',' << means[1] << ',' << departureMax
            << ',' << lambdaMin << ',' << cutoffMax << ',' << reMin << ',' << means[2] << ',' << means[3];
    if (m_diagnostics.is_open())
    {
        for (int count : reasons) m_diagnostics << ',' << count;
        m_diagnostics
            << '\n' << std::flush;
    }
}
} // namespace Nektar
