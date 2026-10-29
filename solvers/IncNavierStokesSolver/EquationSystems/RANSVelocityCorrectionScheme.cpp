///////////////////////////////////////////////////////////////////////////////
// File: RANSVelocityCorrectionScheme.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSVelocityCorrectionScheme.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/ContField.h>
#include <SolverUtils/Forcing/ForcingBody.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>

namespace Nektar
{
using namespace SolverUtils;
using MultiRegions::DirCartesianMap;
std::string RANSVelocityCorrectionScheme::className =
    GetEquationSystemFactory().RegisterCreatorFunction(
        "RANSVelocityCorrectionScheme", RANSVelocityCorrectionScheme::create);
std::string RANSVelocityCorrectionScheme::solverTypeLookupId =
    LibUtilities::SessionReader::RegisterEnumValue(
        "SolverType", "RANSVelocityCorrectionScheme", eRANSVelocityCorrectionScheme);

RANSVelocityCorrectionScheme::RANSVelocityCorrectionScheme(
    const LibUtilities::SessionReaderSharedPtr &session,
    const SpatialDomains::MeshGraphSharedPtr &graph)
    : UnsteadySystem(session, graph), VCSWeakPressure(session, graph)
{
}

NekDouble RANSVelocityCorrectionScheme::Maximum(NekDouble value) const
{
    m_session->GetComm()->AllReduce(value, LibUtilities::ReduceMax);
    return value;
}

void RANSVelocityCorrectionScheme::RequireEverywhere(
    bool valid, const std::string &message) const
{
    int bad = valid ? 0 : 1;
    m_session->GetComm()->AllReduce(bad, LibUtilities::ReduceMax);
    ASSERTL0(!bad, "RANS: " + message);
}

void RANSVelocityCorrectionScheme::v_InitObject(bool declareField)
{
    VCSWeakPressure::v_InitObject(declareField);
    const int dim = m_velocity.size();
    m_k = dim;
    m_epsilon = dim + 1;
    ASSERTL0(((dim == 2 && m_fields[0]->GetExpType() == MultiRegions::e2D) ||
              (dim == 3 && m_fields[0]->GetExpType() == MultiRegions::e3D)) &&
                 m_fields.size() == dim + 3,
             "RANS requires full 2D or 3D CG fields u,v[,w],k,eps,p");
    for (int i = 0; i < dim; ++i)
    {
        ASSERTL0(m_velocity[i] == i, "RANS velocity fields must come first");
    }
    std::string model, diffusion, wallModel;
    m_session->LoadSolverInfo("TurbulenceModel", model, "StandardKEpsilon");
    ASSERTL0(model == "StandardKEpsilon" || model == "CLSHighRe",
             "TurbulenceModel must be StandardKEpsilon or CLSHighRe");
    m_model = model == "CLSHighRe" ? RANSModel::CLSHighRe : RANSModel::StandardKEpsilon;
    ASSERTL0(m_session->GetVariable(m_k) == "k" &&
                 m_session->GetVariable(m_epsilon) ==
                     (m_model == RANSModel::CLSHighRe ? "eps_tilde" : "eps"),
             "RANS fields require k,eps for standard or k,eps_tilde for CLS");
    for (int i = 0; i < m_fields.size(); ++i)
    {
        const auto &field = m_fields[i];
        ASSERTL0(std::dynamic_pointer_cast<MultiRegions::ContField>(field) &&
                     (i == m_fields.size() - 1 || field->GetNcoeffs() == m_fields[0]->GetNcoeffs()) &&
                     field->GetTotPoints() == m_fields[0]->GetTotPoints(),
                 "RANS requires matching u/k/eps expansions and a CG pressure space on the same quadrature grid");
    }
    const int nq = m_fields[0]->GetTotPoints();
    Array<OneD, NekDouble> x(nq, 0.0), y(nq, 0.0), z(nq, 0.0);
    Array<OneD, NekDouble> xf(nq, 0.0), yf(nq, 0.0), zf(nq, 0.0);
    m_fields[0]->GetCoords(x, y, z);
    bool sameGrid = true;
    for (int i = 1; i < m_fields.size(); ++i)
    {
        m_fields[i]->GetCoords(xf, yf, zf);
        for (int q = 0; q < nq; ++q)
            sameGrid &= std::abs(x[q] - xf[q]) < 1e-12 * std::max(1.0, std::abs(x[q])) &&
                        std::abs(y[q] - yf[q]) < 1e-12 * std::max(1.0, std::abs(y[q])) &&
                        std::abs(z[q] - zf[q]) < 1e-12 * std::max(1.0, std::abs(z[q]));
    }
    RequireEverywhere(sameGrid, "all fields must use the same ordered physical quadrature grid");
    ASSERTL0(m_extrapolation->GetSubStepName().empty() && m_cflSafetyFactor == 0 &&
                 m_TimeIncrementFactor == 1,
             "RANS currently requires constant timesteps without substepping");
    for (int i = 0; i < m_nConvectiveFields; ++i)
        ASSERTL0(!m_session->DefinesFunction("DiffusionCoefficient", m_session->GetVariable(i)),
                 "Use RANS transport constants instead of the legacy DiffusionCoefficient function");
    ASSERTL0(m_intScheme->GetName().find("IMEX") == 0 &&
                 m_intScheme->GetVariant().empty() &&
                 (m_intScheme->GetOrder() == 1 || m_intScheme->GetOrder() == 2),
             "RANS supports constant-step IMEX order 1 or 2");
    ASSERTL0(!m_useSpecVanVisc && !m_useHomo1DSpecVanVisc &&
                 !m_useGJPStabilisation && !m_SmoothAdvection && m_flowrate == 0,
             "RANS does not yet support SVV, GJP, smoothing or flowrate control");
    ASSERTL0(!DefinedForcing("MovingReferenceFrame"),
             "RANS requires a stationary Cartesian frame");
    m_session->LoadSolverInfo("RANSResolvedDiffusion", diffusion, "Molecular");
    ASSERTL0(diffusion == "HighRe" || diffusion == "Molecular",
             "RANSResolvedDiffusion must be HighRe or Molecular");
    ASSERTL0(std::isfinite(m_kinvis) && m_kinvis > 0,
             "Kinvis is the positive physical viscosity, also in HighRe mode");
    m_resolvedNu = diffusion == "HighRe" ? 0 : m_kinvis;
    m_session->LoadParameter("RANSCmu", m_parameters.cMu, .09);
    m_session->LoadParameter("RANSCepsilon1", m_parameters.cEpsilon1, 1.44);
    m_session->LoadParameter("RANSCepsilon2", m_parameters.cEpsilon2, 1.92);
    m_session->LoadParameter("RANSSigmaK", m_parameters.sigmaK, 1.0);
    m_session->LoadParameter("RANSSigmaEpsilon", m_parameters.sigmaEpsilon, 1.3);
    m_session->LoadParameter("RANSYapCoefficient", m_parameters.yapCoefficient, .83);
    m_session->LoadParameter("RANSYapLength", m_parameters.yapLength, 2.5);
    m_parameters.Validate();
    m_reference = std::make_shared<const RANSNeutralReference>(
        m_model == RANSModel::CLSHighRe ? RANSNeutralReference::CLS(m_parameters) :
                                        RANSNeutralReference::Standard(m_parameters));
    m_session->LoadParameter("RANSWallOrigin", m_wallOrigin, 0.0);
    ASSERTL0(std::isfinite(m_wallOrigin), "RANSWallOrigin must be finite");
    m_wallDistance = Array<OneD, NekDouble>(nq);
    bool positiveDistance = true;
    for (int q = 0; q < nq; ++q)
    {
        m_wallDistance[q] = (dim == 2 ? y[q] : z[q]) - m_wallOrigin;
        positiveDistance &= m_wallDistance[q] > 0;
    }
    if (m_model == RANSModel::CLSHighRe)
        RequireEverywhere(positiveDistance, "CLS requires positive geometric distance above the flat physical wall");
    m_nonlinearStress = Gradients(dim);
    for (int i = 0; i < dim; ++i)
    {
        m_nonlinearStress[i] = Fields(dim);
        for (int j = 0; j < dim; ++j)
            m_nonlinearStress[i][j] = Array<OneD, NekDouble>(nq, 0.0);
    }
    m_session->LoadParameter("RANSNonlinearTolerance", m_tolerance, 1e-9);
    m_session->LoadParameter("RANSPressureTolerance", m_pressureTolerance, m_tolerance);
    const bool cls = m_model == RANSModel::CLSHighRe;
    m_session->LoadParameter("RANSMaxIterations", m_maxIterations, cls ? 800 : 200);
    m_session->LoadParameter("RANSPressureRelaxation", m_pressureRelaxation, cls ? .25 : .8);
    m_session->LoadParameter("RANSPicardRelaxation", m_relaxation, cls ? .1 : .5);
    m_session->LoadParameter("RANSMinimumK", m_minimumK, 0.0);
    m_session->LoadParameter("RANSMinimumEpsilon", m_minimumEpsilon, 0.0);
    m_session->LoadParameter("RANSZ0", m_z0, .02);
    m_session->LoadParameter("RANSHeightShift", m_heightShift, 0.0);
    ASSERTL0(std::isfinite(m_tolerance) && m_tolerance > 0 &&
                 std::isfinite(m_pressureTolerance) && m_pressureTolerance > 0 &&
                 m_maxIterations > 0 && m_pressureRelaxation > 0 &&
                 m_pressureRelaxation <= 1 && m_relaxation > 0 && m_relaxation <= 1 &&
                 std::isfinite(m_minimumK) && m_minimumK >= 0 &&
                 std::isfinite(m_minimumEpsilon) && m_minimumEpsilon >= 0 &&
                 std::isfinite(m_z0) && m_z0 > 0 &&
                 std::isfinite(m_heightShift) && m_heightShift >= 0,
             "Invalid RANS iteration, positivity or wall parameters");
    m_session->LoadSolverInfo("RANSWallModel", wallModel, "None");
    ASSERTL0(wallModel == "None" || wallModel == "AWFOriginal",
             "RANSWallModel currently supports None or AWFOriginal");
    ASSERTL0(m_model != RANSModel::CLSHighRe || m_heightShift == 0,
             "CLS uses geometric wall distance without a roughness height shift");
    ASSERTL0(wallModel != "AWFOriginal" || m_model == RANSModel::CLSHighRe,
             "The original AWF strip currently requires CLSHighRe");
    bool valid = true;
    for (int b = 0; b < m_pressure->GetBndConditions().size(); ++b)
    {
        // Pressure is a Lagrange multiplier. These markers select the coupled
        // normal-momentum treatment, not a user-prescribed zero pressure flux.
        const auto pbc = m_pressure->GetBndConditions()[b];
        if (pbc->GetBoundaryConditionType() == SpatialDomains::ePeriodic)
            continue;
        valid &= pbc->GetBoundaryConditionType() == SpatialDomains::eNeumann &&
                 pbc->GetUserDefined() == "RANSPressure";
        const auto tag = m_fields[0]->GetBndConditions()[b]->GetUserDefined();
        if (tag == "AWFOriginal")
        {
            valid &= wallModel == "AWFOriginal";
            for (int i = 0; i < dim - 1; ++i)
                valid &= m_fields[i]->GetBndConditions()[b]->GetBoundaryConditionType() == SpatialDomains::eRobin &&
                         m_fields[i]->GetBndConditions()[b]->GetUserDefined() == tag;
            valid &= m_fields[dim - 1]->GetBndConditions()[b]->GetBoundaryConditionType() == SpatialDomains::eDirichlet &&
                     m_fields[m_k]->GetBndConditions()[b]->GetBoundaryConditionType() == SpatialDomains::eRobin &&
                     m_fields[m_k]->GetBndConditions()[b]->GetUserDefined() == tag &&
                     m_fields[m_epsilon]->GetBndConditions()[b]->GetBoundaryConditionType() == SpatialDomains::eDirichlet;
        }
        else if (tag == "RANSLogLaw")
        {
            valid &= m_model == RANSModel::StandardKEpsilon;
            for (int i = 0; i < dim - 1; ++i)
            {
                valid &= m_fields[i]->GetBndConditions()[b]->GetBoundaryConditionType() ==
                             SpatialDomains::eNeumann &&
                         m_fields[i]->GetBndConditions()[b]->GetUserDefined() == tag;
            }
            valid &= m_fields[dim - 1]->GetBndConditions()[b]->GetBoundaryConditionType() == SpatialDomains::eDirichlet &&
                     m_fields[m_k]->GetBndConditions()[b]->GetBoundaryConditionType() == SpatialDomains::eNeumann &&
                     m_fields[m_epsilon]->GetBndConditions()[b]->GetBoundaryConditionType() == SpatialDomains::eDirichlet;
        }
        else
        {
            // Prescribed velocity (including matched inlet/outlet data) or a
            // Cartesian plane with prescribed normal velocity and shear flux.
            Array<OneD, Array<OneD, NekDouble>> normals;
            m_fields[0]->GetBoundaryNormals(b, normals);
            for (int q = 0; q < normals[0].size(); ++q)
                for (int i = 0; i < dim; ++i)
                    if (std::abs(normals[i][q]) > 1e-10)
                        valid &= m_fields[i]->GetBndConditions()[b]->GetBoundaryConditionType() == SpatialDomains::eDirichlet;
        }
    }
    RequireEverywhere(valid, "pressure boundaries need RANSPressure and prescribed normal velocity; log walls need tangential Neumann, normal Dirichlet, k Neumann, eps Dirichlet");
    if (wallModel == "AWFOriginal")
    {
        for (const auto &force : m_forcing)
            ASSERTL0(std::dynamic_pointer_cast<ForcingBody>(force),
                "The initial AWF adapter supports Body/Field forcing only, since its current forcing is evaluated within the implicit iteration");
        m_wallStrip = std::make_shared<RANSAnalyticalWallStrip>(
            m_session, m_fields, dim, m_wallOrigin);
    }
    // Pressure's Poisson map has an artificial pinned vertex. Mass inverses
    // used by the preconditioner must span the FULL pressure space, including
    // that vertex: a mass operator is nonsingular without any pin.
    m_pressureMass = MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(
        m_session, m_fields[0]->GetGraph(), "p", true, false);
}

void RANSVelocityCorrectionScheme::CheckState(const Fields &state) const
{
    bool valid = true;
    double minK = std::numeric_limits<double>::max(), minE = minK;
    for (int f = 0; f < state.size(); ++f)
        for (int q = 0; q < state[f].size(); ++q)
            valid &= std::isfinite(state[f][q]);
    for (int q = 0; q < state[m_k].size(); ++q)
    {
        valid &= state[m_k][q] > m_minimumK &&
                 state[m_epsilon][q] > m_minimumEpsilon;
        minK = std::min(minK, state[m_k][q]);
        minE = std::min(minE, state[m_epsilon][q]);
    }
    m_session->GetComm()->AllReduce(minK, LibUtilities::ReduceMin);
    m_session->GetComm()->AllReduce(minE, LibUtilities::ReduceMin);
    RequireEverywhere(valid, "nonfinite state or k/eps positivity guard activated; min k=" +
        std::to_string(minK) + ", min eps=" + std::to_string(minE) + "; no floor is applied");
}

void RANSVelocityCorrectionScheme::v_DoInitialise(bool dump)
{
    VCSWeakPressure::v_DoInitialise(dump);
    Fields initial(m_nConvectiveFields);
    for (int i = 0; i < initial.size(); ++i)
        initial[i] = m_fields[i]->UpdatePhys();
    CheckState(initial);
    m_fieldMetaDataMap["RANSModel"] = m_model == RANSModel::CLSHighRe ? "CLSHighRe" : "StandardKEpsilon";
    m_fieldMetaDataMap["RANSDissipation"] = "high-Re: eps_total=eps_transport; molecular mode adds diffusion only";
    m_fieldMetaDataMap["RANSPhysicalFrame"] = m_velocity.size() == 2 ? "mesh(x,y)=physical(x,z)" : "physical(x,y,z)";
    m_fieldMetaDataMap["RANSWallModel"] = m_wallStrip ? "AWFOriginal conservative strip" : "None or standard log law";
    m_fieldMetaDataMap["RANSPressure"] = "pi=p/rho+2k/3";
    m_fieldMetaDataMap["RANSResolvedNu"] = std::to_string(m_resolvedNu);
    if (m_session->GetComm()->GetRank() == 0)
    {
        m_diagnostics.open(m_session->GetSessionName() + ".rans.csv");
        ASSERTL0(m_diagnostics.good(), "Cannot open RANS diagnostics");
        m_diagnostics << "time,iterations,change,pressure_correction,divergence_l2,weak_divergence_l2,min_k,min_eps,min_nut,max_nut,live_global_systems,guard_activations\n" << std::setprecision(17);
    }
}

void RANSVelocityCorrectionScheme::v_GenerateSummary(SummaryList &summary)
{
    VCSWeakPressure::v_GenerateSummary(summary);
    AddSummaryItem(summary, "RANS model", m_model == RANSModel::CLSHighRe ?
        "High-Re CLS, full cubic stress and algebraic Yap" : "Standard k-epsilon, full symmetric stress");
    AddSummaryItem(summary, "RANS resolved molecular viscosity", m_resolvedNu);
    AddSummaryItem(summary, "RANS neutral kappa", m_reference->Kappa());
    AddSummaryItem(summary, "RANS implicit stage", "nonlinear Picard + pressure defect iteration");
    AddSummaryItem(summary, "RANS pressure correction tolerance", m_pressureTolerance);
}

void RANSVelocityCorrectionScheme::v_GetPressure(
    const Array<OneD, const Array<OneD, NekDouble>> &fields,
    Array<OneD, NekDouble> &pressure)
{
    pressure = Array<OneD, NekDouble>(m_pressure->GetTotPoints());
    m_pressure->BwdTrans(m_pressure->GetCoeffs(), pressure);
    Vmath::Svtvp(pressure.size(), -2.0 / 3.0, fields[m_k], 1,
                 pressure, 1, pressure, 1);
}

void RANSVelocityCorrectionScheme::v_ExtraFldOutput(
    std::vector<Array<OneD, NekDouble>> &coefficients,
    std::vector<std::string> &variables)
{
    if (m_model != RANSModel::CLSHighRe) return;
    const int dim = m_velocity.size(), nq = m_fields[0]->GetTotPoints();
    Gradients gradient(dim);
    Fields state(m_nConvectiveFields);
    for (int i = 0; i < state.size(); ++i)
    {
        state[i] = Array<OneD, NekDouble>(nq);
        m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), state[i]);
    }
    for (int i = 0; i < dim; ++i)
    {
        gradient[i] = Fields(dim);
        for (int j = 0; j < dim; ++j)
        {
            gradient[i][j] = Array<OneD, NekDouble>(nq);
            m_fields[i]->PhysDeriv(DirCartesianMap[j], state[i], gradient[i][j]);
        }
    }
    const std::vector<std::string> names = {"rans_nut", "rans_P", "rans_Yap", "rans_eps_total",
        "rans_R_xx", "rans_R_yy", "rans_R_zz", "rans_R_xy", "rans_R_xz", "rans_R_yz",
        "rans_quadratic_norm", "rans_cubic_norm"};
    Fields values(names.size());
    for (auto &value : values) value = Array<OneD, NekDouble>(nq);
    bool valid = true;
    try
    {
        for (int q = 0; q < nq; ++q)
        {
            RANSTensor g{};
            const auto axis = [dim](int i) { return dim == 2 && i == 1 ? 2 : i; };
            for (int i = 0; i < dim; ++i)
                for (int j = 0; j < dim; ++j) g[axis(i)][axis(j)] = gradient[i][j][q];
            const auto v = EvaluateRANSClosure(m_model, state[m_k][q], state[m_epsilon][q],
                g, m_wallDistance[q], m_parameters, m_resolvedNu);
            values[0][q] = v.eddyViscosity; values[1][q] = v.production;
            values[2][q] = v.yap; values[3][q] = state[m_epsilon][q];
            const std::array<std::array<int, 2>, 6> entries{{{0,0},{1,1},{2,2},{0,1},{0,2},{1,2}}};
            for (int c = 0; c < entries.size(); ++c)
                values[4+c][q] = v.reynoldsStress[entries[c][0]][entries[c][1]];
            double quadratic = 0, cubic = 0;
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                {
                    quadratic += v.quadraticStress[i][j] * v.quadraticStress[i][j];
                    cubic += v.cubicStress[i][j] * v.cubicStress[i][j];
                }
            values[10][q] = std::sqrt(quadratic); values[11][q] = std::sqrt(cubic);
        }
    }
    catch (const std::exception &) { valid = false; }
    RequireEverywhere(valid, "cannot output an inadmissible CLS state");
    for (int i = 0; i < names.size(); ++i)
    {
        Array<OneD, NekDouble> data(m_fields[0]->GetNcoeffs());
        m_fields[0]->FwdTransLocalElmt(values[i], data);
        coefficients.push_back(data); variables.push_back(names[i]);
    }
}

void RANSVelocityCorrectionScheme::v_EvaluateAdvection_SetPressureBCs(
    const Array<OneD, const Array<OneD, NekDouble>> &input,
    Fields &output, NekDouble time)
{
    EvaluateAdvectionTerms(input, output, time);
    for (auto &force : m_forcing)
        force->Apply(m_fields, input, output, time);
    // No legacy curl-curl pressure extrapolation. The implicit stage enforces
    // continuity with the actual variable symmetric stress and stage BCs.
}

void RANSVelocityCorrectionScheme::UpdateLogWall(const Fields &state)
{
    const int dim = m_velocity.size();
    bool valid = true;
    for (int b = 0; b < m_fields[0]->GetBndConditions().size(); ++b)
    {
        if (m_fields[0]->GetBndConditions()[b]->GetUserDefined() != "RANSLogLaw")
            continue;
        const auto exp = m_fields[0]->GetBndCondExpansions()[b];
        const int nq = exp->GetTotPoints();
        Fields values(dim), normals;
        Array<OneD, NekDouble> k, x(nq, 0.0), y(nq, 0.0), z(nq, 0.0), eps(nq);
        exp->GetCoords(x, y, z);
        m_fields[0]->GetBoundaryNormals(b, normals);
        m_fields[m_k]->ExtractPhysToBnd(b, state[m_k], k);
        for (int i = 0; i < dim; ++i)
            m_fields[i]->ExtractPhysToBnd(b, state[i], values[i]);
        for (int q = 0; q < nq; ++q)
        {
            const double height = dim == 2 ? y[q] : z[q];
            valid &= std::abs(normals[dim - 1][q] + 1) < 1e-10 &&
                     std::abs(values[dim - 1][q]) < 1e-10 &&
                     height + m_heightShift > m_z0 && k[q] > m_minimumK;
        }
        // Check collectively after the local loop; calls below need valid
        // geometry, checked without throwing on a single partition first.
        if (!valid) continue;
        for (int q = 0; q < nq; ++q)
        {
            const double height = dim == 2 ? y[q] : z[q];
            const auto tau = m_reference->StandardLogWallTraction(
                k[q], {values[0][q], dim == 3 ? values[1][q] : 0.0},
                height, m_z0, m_heightShift);
            for (int i = 0; i < dim - 1; ++i) values[i][q] = -tau[i];
            const double uk = std::pow(m_parameters.cMu, .25) * std::sqrt(k[q]);
            eps[q] = uk * uk * uk / (m_reference->Kappa() * (height + m_heightShift));
        }
        for (int i = 0; i < dim - 1; ++i)
        {
            auto boundary = m_fields[i]->GetBndCondExpansions()[b];
            boundary->IProductWRTBase(values[i], boundary->UpdateCoeffs());
        }
        auto boundary = m_fields[m_epsilon]->GetBndCondExpansions()[b];
        boundary->FwdTransBndConstrained(eps, boundary->UpdateCoeffs());
    }
    RequireEverywhere(valid, "RANSLogLaw requires a flat lower wall, zero normal velocity, positive k and height+shift > z0");
}

void RANSVelocityCorrectionScheme::Closure(
    const Fields &state, Gradients &gradient, Fields &diffusion,
    Fields &gain, Fields &loss)
{
    const int nq = m_fields[0]->GetTotPoints(), dim = m_velocity.size();
    CheckState(state);
    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
            m_fields[i]->PhysDeriv(DirCartesianMap[j], state[i], gradient[i][j]);
    bool valid = true;
    try
    {
      for (int q = 0; q < nq; ++q)
      {
        RANSTensor g{};
        const auto axis = [dim](int i) { return dim == 2 && i == 1 ? 2 : i; };
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j) g[axis(i)][axis(j)] = gradient[i][j][q];
        const auto v = EvaluateRANSClosure(m_model,
            state[m_k][q], state[m_epsilon][q], g, m_wallDistance[q], m_parameters, m_resolvedNu);
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j)
                m_nonlinearStress[i][j][q] = v.quadraticStress[axis(i)][axis(j)] + v.cubicStress[axis(i)][axis(j)];
        for (int i = 0; i < dim; ++i) diffusion[i][q] = v.momentumDiffusion;
        diffusion[m_k][q] = v.kDiffusion;
        diffusion[m_epsilon][q] = v.epsilonDiffusion;
        gain[m_k][q] = v.production;
        loss[m_k][q] = state[m_epsilon][q] / state[m_k][q];
        gain[m_epsilon][q] = v.epsilonProduction + v.yap;
        loss[m_epsilon][q] = m_parameters.cEpsilon2 * loss[m_k][q];
      }
    }
    catch (const std::exception &) { valid = false; }
    RequireEverywhere(valid, "point closure rejected a nonfinite gradient or arithmetic range failure");
}

void RANSVelocityCorrectionScheme::SolveField(
    int field, const Array<OneD, NekDouble> &rhs,
    const Array<OneD, NekDouble> &diffusion,
    const Array<OneD, NekDouble> &mass, Array<OneD, NekDouble> &output)
{
    StdRegions::VarCoeffMap coefficients;
    for (auto type : {StdRegions::eVarCoeffD00, StdRegions::eVarCoeffD11})
        coefficients[type] = diffusion;
    if (m_velocity.size() == 3) coefficients[StdRegions::eVarCoeffD22] = diffusion;
    coefficients[StdRegions::eVarCoeffMass] = mass;
    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorLambda] = 1.0;
    const auto key = m_fields[field]->HelmSolve(
        rhs, m_fields[field]->UpdateCoeffs(), factors, coefficients,
        StdRegions::NullVarFactorsMap, NullNekDouble1DArray, false);
    m_fields[field]->BwdTrans(m_fields[field]->GetCoeffs(), output);
    // Every changing-coefficient solve releases its global, static-condensed
    // and local blocks. No coefficient history accumulates in matrix caches.
    m_fields[field]->UnsetGlobalLinSys(key, true);
}

NekDouble RANSVelocityCorrectionScheme::CorrectPressure(
    const Fields &velocity, const Array<OneD, NekDouble> &mu,
    Array<OneD, NekDouble> &pressure, NekDouble alpha, NekDouble &divergence,
    NekDouble &weakDivergence)
{
    const int nq = m_fields[0]->GetTotPoints(), nc = m_pressure->GetNcoeffs();
    const int dim = m_velocity.size();
    Array<OneD, NekDouble> rhs(nc), weak(nc), increment(nq), div(nq, 0.0), tmp(nq), projection(nc, 0.0);
    m_pressure->IProductWRTDerivBase(velocity, rhs);
    Vmath::Smul(nc, -1.0, rhs, 1, weak, 1);
    Vmath::Smul(nc, -alpha, rhs, 1, rhs, 1);
    // Boundary normal velocity completes the weak divergence. It is evaluated
    // from the current CG trace, whose normal component is prescribed. This
    // replaces the constant-viscosity high-order boundary extrapolation.
    for (int b = 0; b < m_pressure->GetBndConditions().size(); ++b)
    {
        if (m_pressure->GetBndConditions()[b]->GetBoundaryConditionType() == SpatialDomains::ePeriodic)
            continue;
        Fields trace(dim);
        for (int i = 0; i < dim; ++i)
            m_fields[i]->ExtractPhysToBnd(b, velocity[i], trace[i]);
        auto boundary = m_pressure->GetBndCondExpansions()[b];
        boundary->NormVectorIProductWRTBase(trace, boundary->UpdateCoeffs());
    }
    m_pressure->ImposeNeumannConditions(weak);
    for (int b = 0; b < m_pressure->GetBndCondExpansions().size(); ++b)
    {
        if (m_pressure->GetBndConditions()[b]->GetBoundaryConditionType() == SpatialDomains::ePeriodic)
            continue;
        auto boundary = m_pressure->GetBndCondExpansions()[b];
        Vmath::Smul(boundary->GetNcoeffs(), -alpha, boundary->GetCoeffs(), 1,
                     boundary->UpdateCoeffs(), 1);
    }
    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorLambda] = 0;
    // The pure-Neumann assembly map pins a vertex using the supplied initial
    // coefficients. An increment has a homogeneous pin, not the old pi value.
    Vmath::Zero(nc, m_pressure->UpdateCoeffs(), 1);
    m_pressure->HelmSolve(rhs, m_pressure->UpdateCoeffs(), factors,
        StdRegions::NullVarCoeffMap, StdRegions::NullVarFactorsMap,
        NullNekDouble1DArray, false);
    m_pressure->BwdTrans(m_pressure->GetCoeffs(), increment);
    for (int i = 0; i < dim; ++i)
    {
        m_fields[i]->PhysDeriv(DirCartesianMap[i], velocity[i], tmp);
        Vmath::Vadd(nq, tmp, 1, div, 1, div, 1);
    }
    divergence = m_fields[0]->L2(div);
    // Cahouet--Chabard-type pressure preconditioning: the Poisson part above
    // handles inertia, and -mu div(u) handles the viscous limit. Only the
    // preconditioner is approximate; convergence tests the coupled equations.
    // Apply the viscous preconditioner to the SAME assembled weak divergence,
    // rather than a differently weighted strong divergence. This distinction
    // matters for variable mu and a mixed velocity/pressure approximation.
    m_pressureMass->MultiplyByInvMassMatrix(weak, projection);
    m_pressureMass->BwdTrans(projection, tmp);
    weakDivergence = m_pressureMass->L2(tmp);
    Vmath::Vmul(nq, mu, 1, tmp, 1, tmp, 1);
    m_pressureMass->FwdTrans(tmp, projection);
    m_pressureMass->BwdTrans(projection, tmp);
    Vmath::Vsub(nq, increment, 1, tmp, 1, increment, 1);
    Array<OneD, NekDouble> ones(nq, 1.0);
    const double mean = m_pressure->Integral(increment) / m_pressure->Integral(ones);
    NekDouble change = 0;
    for (int q = 0; q < nq; ++q)
    {
        increment[q] -= mean;
        pressure[q] += m_pressureRelaxation * increment[q];
        change = std::max(change, std::abs(increment[q]));
    }
    return Maximum(change);
}

void RANSVelocityCorrectionScheme::v_SolveUnsteadyStokesSystem(
    const Array<OneD, const Array<OneD, NekDouble>> &input,
    Fields &output, NekDouble time, NekDouble lambda)
{
    ASSERTL0(lambda > 0, "RANS stage requires positive implicit timestep");
    const int nf = m_nConvectiveFields, nq = m_fields[0]->GetTotPoints();
    const int nc = m_fields[0]->GetNcoeffs(), dim = m_velocity.size();
    const double alpha = 1.0 / lambda;
    Fields state(nf), next(nf), diffusion(nf), gain(nf), loss(nf);
    Gradients gradient(dim);
    for (int i = 0; i < nf; ++i)
    {
        state[i] = Array<OneD, NekDouble>(nq);
        next[i] = Array<OneD, NekDouble>(nq);
        diffusion[i] = Array<OneD, NekDouble>(nq);
        gain[i] = Array<OneD, NekDouble>(nq, 0.0);
        loss[i] = Array<OneD, NekDouble>(nq, 0.0);
        Vmath::Vcopy(nq, m_fields[i]->GetPhys(), 1, state[i], 1);
    }
    for (int i = 0; i < dim; ++i)
    {
        gradient[i] = Fields(dim);
        for (int j = 0; j < dim; ++j) gradient[i][j] = Array<OneD, NekDouble>(nq);
    }
    Array<OneD, NekDouble> pressure(nq), rhs(nc), cross(nc), tmp(nq), mass(nq), gradP(nq);
    m_pressure->BwdTrans(m_pressure->GetCoeffs(), pressure);
    Array<OneD, NekDouble> ones(nq, 1.0);
    const double pressureMean = m_pressure->Integral(pressure) / m_pressure->Integral(ones);
    for (int q = 0; q < nq; ++q) pressure[q] -= pressureMean;
    SetBoundaryConditions(time);
    double change = 0, correction = 0, div = 0, weakDiv = 0;
    int iteration = 0;
    for (; iteration < m_maxIterations; ++iteration)
    {
        Closure(state, gradient, diffusion, gain, loss);
        UpdateLogWall(state);
        if (m_wallStrip)
        {
            Fields advection(nf), body(nf), pressureForce(dim), acceleration(dim);
            for (int i = 0; i < nf; ++i)
            {
                advection[i] = Array<OneD, NekDouble>(nq, 0.0);
                body[i] = Array<OneD, NekDouble>(nq, 0.0);
            }
            EvaluateAdvectionTerms(state, advection, time);
            for (auto &force : m_forcing) force->Apply(m_fields, state, body, time);
            Array<OneD, NekDouble> mechanical(nq);
            Vmath::Svtvp(nq, -2.0/3.0, state[m_k], 1, pressure, 1, mechanical, 1);
            for (int i = 0; i < dim; ++i)
            {
                pressureForce[i] = Array<OneD, NekDouble>(nq);
                acceleration[i] = Array<OneD, NekDouble>(nq);
                m_fields[0]->PhysDeriv(DirCartesianMap[i], mechanical, pressureForce[i]);
                Vmath::Smul(nq, -1.0, advection[i], 1, acceleration[i], 1);
                Vmath::Vsub(nq, acceleration[i], 1, body[i], 1, acceleration[i], 1);
            }
            Fields predictor(nf);
            for (int i = 0; i < nf; ++i)
            {
                predictor[i] = Array<OneD, NekDouble>(nq);
                Vmath::Vcopy(nq, input[i], 1, predictor[i], 1);
            }
            m_wallStrip->Update(state, predictor, gradient, pressureForce, acceleration, alpha);
        }
        for (int i = 0; i < nf; ++i)
        {
            for (int q = 0; q < nq; ++q)
            {
                tmp[q] = -alpha * input[i][q] - gain[i][q];
                mass[q] = alpha + loss[i][q];
            }
            if (i < dim)
            {
                m_pressure->PhysDeriv(DirCartesianMap[i], pressure, gradP);
                Vmath::Vadd(nq, gradP, 1, tmp, 1, tmp, 1);
            }
            m_fields[i]->IProductWRTBase(tmp, rhs);
            if (i < dim)
            {
                Fields flux(dim);
                for (int j = 0; j < dim; ++j)
                {
                    flux[j] = Array<OneD, NekDouble>(nq);
                    Vmath::Vmul(nq, diffusion[i], 1, gradient[j][i], 1, flux[j], 1);
                    Vmath::Vsub(nq, flux[j], 1, m_nonlinearStress[i][j], 1, flux[j], 1);
                }
                m_fields[i]->IProductWRTDerivBase(flux, cross);
                Vmath::Vadd(nc, cross, 1, rhs, 1, rhs, 1);
            }
            SolveField(i, rhs, diffusion[i], mass, next[i]);
            if (m_session->DefinesParameter("RANSDebug") && m_session->GetComm()->GetRank() == 0)
            {
                std::cout << "RANS iteration " << iteration << " field " << i
                          << " input " << *std::min_element(input[i].begin(), input[i].end())
                          << " gain " << *std::min_element(gain[i].begin(), gain[i].end())
                          << " result [" << *std::min_element(next[i].begin(), next[i].end())
                          << ',' << *std::max_element(next[i].begin(), next[i].end()) << "]" << std::endl;
            }
        }
        CheckState(next);
        Fields velocity(dim);
        for (int i = 0; i < dim; ++i) velocity[i] = next[i];
        correction = CorrectPressure(velocity, diffusion[0], pressure, alpha, div, weakDiv);
        if (m_session->DefinesParameter("RANSDebug") && m_session->GetComm()->GetRank() == 0)
            std::cout << "RANS pressure " << iteration << " correction " << correction
                      << " divergence " << div << " weak " << weakDiv << " range ["
                      << *std::min_element(pressure.begin(), pressure.end()) << ','
                      << *std::max_element(pressure.begin(), pressure.end()) << ']' << std::endl;
        change = 0;
        for (int i = 0; i < nf; ++i)
        {
            double scale = 1e-12, difference = 0;
            for (int q = 0; q < nq; ++q)
            {
                scale = std::max(scale, std::abs(next[i][q]));
                difference = std::max(difference, std::abs(next[i][q] - state[i][q]));
            }
            scale = Maximum(scale);
            change = std::max(change, Maximum(difference) / std::max(1.0, scale));
            for (int q = 0; q < nq; ++q)
                state[i][q] += m_relaxation * (next[i][q] - state[i][q]);
        }
        double pressureScale = 1;
        for (double p : pressure) pressureScale = std::max(pressureScale, std::abs(p));
        if (change < m_tolerance && correction / Maximum(pressureScale) < m_pressureTolerance && weakDiv < m_tolerance)
            break;
    }
    RequireEverywhere(iteration < m_maxIterations, "coupled implicit stage failed to converge");
    // Return the actual Helmholtz solutions. They have the prescribed traces
    // and correspond to the last fixed-point residual, with no profile reset.
    for (int i = 0; i < nf; ++i)
        Vmath::Vcopy(nq, next[i], 1, output[i], 1);
    if (m_wallStrip)
        m_wallStrip->WriteDiagnostics(time, next, gain, loss, pressure);
    // pi is already continuous (a sum of CG increments). A global mass solve
    // on the singular pressure map would impose its pinned vertex on pi.
    m_pressure->FwdTransLocalElmt(pressure, m_pressure->UpdateCoeffs());
    m_pressure->AvgAssemble();
    m_pressure->GlobalToLocal();
    m_pressure->BwdTrans(m_pressure->GetCoeffs(), m_pressure->UpdatePhys());
    double minK = std::numeric_limits<double>::max(), minE = minK, minNu = minK, maxNu = 0;
    int live = 0;
    for (int i = 0; i < nf; ++i) live += m_fields[i]->GetPoolCount("GlobalLinSys");
    for (int q = 0; q < nq; ++q)
    {
        minK = std::min(minK, next[m_k][q]); minE = std::min(minE, next[m_epsilon][q]);
        minNu = std::min(minNu, diffusion[0][q] - m_resolvedNu);
        maxNu = std::max(maxNu, diffusion[0][q] - m_resolvedNu);
    }
    m_session->GetComm()->AllReduce(minK, LibUtilities::ReduceMin);
    m_session->GetComm()->AllReduce(minE, LibUtilities::ReduceMin);
    m_session->GetComm()->AllReduce(minNu, LibUtilities::ReduceMin);
    maxNu = Maximum(maxNu);
    m_session->GetComm()->AllReduce(live, LibUtilities::ReduceMax);
    if (m_diagnostics.is_open())
        m_diagnostics << time << ',' << iteration + 1 << ',' << change << ',' << correction
                      << ',' << div << ',' << weakDiv << ',' << minK << ',' << minE << ',' << minNu << ','
                      << maxNu << ',' << live << ",0\n" << std::flush;
}
} // namespace Nektar
