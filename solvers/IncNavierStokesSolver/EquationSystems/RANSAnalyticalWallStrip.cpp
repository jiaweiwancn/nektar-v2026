///////////////////////////////////////////////////////////////////////////////
// File: RANSAnalyticalWallStrip.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSAnalyticalWallStrip.h>
#include <LibUtilities/Communication/Comm.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace Nektar
{
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
    int dimension, NekDouble wallOrigin)
    : m_session(session), m_fields(fields), m_dim(dimension), m_localPoints(0),
      m_globalPoints(0), m_globalOffset(0), m_remoteSamples(0)
{
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
        if (m_fields[0]->GetBndConditions()[b]->GetUserDefined() == "AWFOriginal")
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
    Require(valid, "AWFOriginal must occupy the flat lower interface at wallOrigin+H");
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
    Require(m_globalPoints > 0, "AWFOriginal requires a tagged wall region");
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
    if (comm->GetRank() == 0)
    {
        m_diagnostics.open(m_session->GetSessionName() + ".awf.csv");
        ASSERTL0(m_diagnostics.good(), "Cannot open AWF strip diagnostics");
        m_diagnostics << std::setprecision(17)
            << "time,area,traction_x,traction_y,mean_production,mean_dissipation,min_kP,min_epsP,epsilonH_error,mass_error,normal_velocity_max,momentum_budget_x,momentum_budget_y,k_budget,budget_available,global_wall_points,remote_N_samples,chi,guard_activations,mean_uH,mean_uN,mean_uP,mean_kP,mean_forcing_x,mean_epsP,mean_epsH,wall_ranks,remote_N_without_wall_samples,epsilonH_projection_error\n";
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

void RANSAnalyticalWallStrip::Update(const Fields &state, const Fields &predictor,
    const Gradients &gradient, const Fields &pressureForce,
    const Fields &convectionMinusBodyForce, NekDouble alpha)
{
    m_volumePredictor = predictor;
    const auto h = Trace(state), n = SampleN(state);
    const auto bh = Trace(predictor), bn = SampleN(predictor);
    const auto ph = Trace(pressureForce), pn = SampleN(pressureForce);
    const auto ah = Trace(convectionMinusBodyForce), an = SampleN(convectionMinusBodyForce);
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
            AnalyticalWallFunction wall(m_awfGeometry, m_kP[q], m_parameters);
            AWFVector velocity{}, forcing{};
            for (int i = 0; i < m_dim - 1; ++i)
            {
                velocity[i] = h[i][q];
                m_pressureP[i][q] = m_geometry.EliminateVelocityP(ph[i][q], pn[i][q]);
                forcing[i] = m_pressureP[i][q] + m_geometry.EliminateVelocityP(ah[i][q], an[i][q]);
            }
            const auto values = wall.Evaluate(velocity, forcing);
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
    Require(valid, "nonpositive eliminated KP or inadmissible/nonfinite original AWF state; no floor is applied");
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
            bc->GetUserDefined() == "AWFOriginal") continue;
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
    if (m_diagnostics.is_open())
        m_diagnostics << time << ',' << area << ',' << tx << ',' << ty << ',' << production << ',' << dissipation
            << ',' << minK << ',' << minE << ',' << epsError << ',' << massError << ',' << normalMax
            << ',' << budget[0] << ',' << budget[1] << ',' << budget[2] << ',' << allAvailable
            << ',' << m_globalPoints << ',' << m_remoteSamples << ",0,0"
            << ',' << uh << ',' << un << ',' << up << ',' << kp << ',' << forcing << ',' << ep << ',' << eh
            << ',' << m_wallRanks << ',' << m_remoteWithoutWallSamples
            << ',' << projectionError
            << '\n' << std::flush;
}
} // namespace Nektar
