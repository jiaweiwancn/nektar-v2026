///////////////////////////////////////////////////////////////////////////////
// File: RANSVelocityCorrectionScheme.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSVelocityCorrectionScheme.h>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/LinearAlgebra/NekLinSysIter.h>
#include <MultiRegions/ContField.h>
#include <SolverUtils/Forcing/ForcingBody.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <filesystem>
#include <sstream>

namespace Nektar
{
using namespace SolverUtils;
using MultiRegions::DirCartesianMap;
namespace
{
struct RANSLinearAction
{
    LibUtilities::NekSysOperators::FunctorType1 action;
    void Apply(LibUtilities::NekSysOperators::InArrayType &input,
               LibUtilities::NekSysOperators::OutArrayType &output,const bool &flag)
    { action(input,output,flag); }
};
struct RANSFlagGuard
{
    bool &flag,old;
    RANSFlagGuard(bool &value,bool active) : flag(value),old(value) { flag = active; }
    ~RANSFlagGuard() { flag = old; }
};
struct RANSTrialGuard
{
    bool &flag;
    RANSWallBoundary *wall;
    RANSTrialGuard(bool &value,RANSWallBoundary *boundary,bool trial) : flag(value),wall(boundary)
    { flag = trial; if (wall) wall->SetTrialMode(trial); }
    ~RANSTrialGuard() { flag = false; if (wall) wall->SetTrialMode(false); }
};
// A positive outlet impedance belongs only to the approximate pressure
// inverse. It retains the constant pressure mode needed to correct net outlet
// flux. The actual multiplier and continuity have no added boundary equation.
class RANSOutletPressureInverse : public MultiRegions::ContField
{
public:
    RANSOutletPressureInverse(const LibUtilities::SessionReaderSharedPtr &session,
        const SpatialDomains::MeshGraphSharedPtr &graph, const std::vector<int> &outlets,
        NekDouble impedance)
        : ContField(session,graph,"p",true,false)
    {
        std::ostringstream coefficient; coefficient << std::setprecision(17) << impedance;
        for (int b : outlets)
        {
            m_bndConditions[b] = std::make_shared<SpatialDomains::RobinBoundaryCondition>(session,"0",coefficient.str());
            Vmath::Zero(m_bndCondExpansions[b]->GetNcoeffs(),m_bndCondExpansions[b]->UpdateCoeffs(),1);
        }
        m_locToGloMap = MemoryManager<MultiRegions::AssemblyMapCG>::AllocateSharedPtr(
            m_session,m_ncoeffs,*this,m_bndCondExpansions,m_bndConditions,false,"p",
            m_periodicVerts,m_periodicEdges,m_periodicFaces);
    }
};
}
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
    if (bad && m_nonlinearTrial) throw std::runtime_error("RANS nonlinear trial: "+message);
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
    m_wallGeometry = ReadRANSWallGeometry(m_session,m_wallOrigin);
    m_terrainWall = m_wallGeometry.terrain != RANSTerrain::Flat;
    std::string advectionTreatment;
    m_session->LoadSolverInfo("RANSAdvectionTreatment",advectionTreatment,"Explicit");
    ASSERTL0(advectionTreatment == "Explicit" || advectionTreatment == "Implicit",
             "RANSAdvectionTreatment must be Explicit or Implicit");
    m_implicitAdvection = advectionTreatment == "Implicit";
    m_session->LoadParameter("RANSGradDiv", m_gradDiv, 0.0);
    ASSERTL0(std::isfinite(m_gradDiv) && m_gradDiv >= 0 &&
             (m_gradDiv == 0 || m_terrainWall),
             "RANSGradDiv must be nonnegative and uses the terrain vector solver");
    m_session->LoadParameter("RANSStartupDuration", m_startupDuration, 0.0);
    ASSERTL0(std::isfinite(m_startupDuration) && m_startupDuration >= 0,
             "RANSStartupDuration must be finite and nonnegative");
    ASSERTL0(m_startupDuration == 0 ||
        (m_wallGeometry.terrain == RANSTerrain::Plane &&
         m_session->DefinesSolverInfo("RANSExperimentalSignature") && m_forcing.empty()),
        "startup continuation is restricted to an experimental plane without external forcing");
    m_wallDistance = Array<OneD, NekDouble>(nq);
    bool positiveDistance = true;
    for (int q = 0; q < nq; ++q)
    {
        m_wallDistance[q] = m_terrainWall ?
            m_wallGeometry.Closest(x[q],dim == 2 ? y[q] : z[q]).distance :
            (dim == 2 ? y[q] : z[q]) - m_wallOrigin;
        positiveDistance &= m_wallDistance[q] > 0;
    }
    if (m_model == RANSModel::CLSHighRe)
        RequireEverywhere(positiveDistance, "CLS requires positive geometric distance above the physical wall");
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
    m_session->LoadParameter("RANSMaxIterations", m_maxIterations, m_terrainWall ? 4000 : (cls ? 800 : 200));
    m_session->LoadParameter("RANSPressureRelaxation", m_pressureRelaxation, m_terrainWall ? .1 : (cls ? .25 : .8));
    m_session->LoadParameter("RANSPicardRelaxation", m_relaxation, m_terrainWall ? .1 : (cls ? .1 : .5));
    std::string nonlinearSolver;
    m_session->LoadSolverInfo("RANSNonlinearSolver",nonlinearSolver,m_terrainWall ? "NewtonKrylov" : "Picard");
    ASSERTL0(nonlinearSolver == "Picard" || nonlinearSolver == "NewtonKrylov",
             "RANSNonlinearSolver must be Picard or NewtonKrylov");
    m_useNewton = nonlinearSolver == "NewtonKrylov";
    std::string newtonMap;
    m_session->LoadSolverInfo("RANSNewtonMap",newtonMap,"CurrentSolve");
    ASSERTL0(newtonMap == "CurrentSolve" || newtonMap == "PreconditionedResidual",
             "RANSNewtonMap must be CurrentSolve or PreconditionedResidual");
    m_preconditionedNewton = newtonMap == "PreconditionedResidual";
    ASSERTL0(!m_preconditionedNewton || (m_useNewton && m_terrainWall),
             "PreconditionedResidual requires terrain NewtonKrylov");
    m_session->LoadParameter("RANSNewtonDifferenceStep",m_newtonDifference,cls ? 1e-6 : 1e-7);
    ASSERTL0(std::isfinite(m_newtonDifference) && m_newtonDifference > 0 && m_newtonDifference < .01,
             "RANSNewtonDifferenceStep must be finite and in (0,0.01)");
    m_session->LoadParameter("RANSNewtonKrylovDimension",m_newtonKrylovDimension,240);
    ASSERTL0(m_newtonKrylovDimension >= 20 && m_newtonKrylovDimension <= 512,
             "RANSNewtonKrylovDimension must be 20..512");
    ASSERTL0(!m_useNewton || m_terrainWall,"NewtonKrylov is supported on the terrain path");
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
    ASSERTL0(m_terrainWall || wallModel != "None" || (!m_session->DefinesFunction("RANSWallProperties","z0") && !m_session->DefinesFunction("RANSWallProperties","roughness")),
             "Select RANSWallGeometry=Plane for spatial log-wall properties");
    ASSERTL0(wallModel == "None" || wallModel == "AWFOriginal" || wallModel == "AWFEquilibrium",
             "RANSWallModel supports None, AWFOriginal or AWFEquilibrium");
    ASSERTL0(m_model != RANSModel::CLSHighRe || m_heightShift == 0,
             "CLS uses geometric wall distance without a roughness height shift");
    ASSERTL0(wallModel == "None" || m_heightShift == 0,
             "The analytical wall strip uses geometric height without RANSHeightShift");
    bool valid = true, terrainBoundary = false;
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
        if (tag == "RANSOutflow")
        {
            m_outflows.push_back(b);
            Fields normals;
            m_fields[0]->GetBoundaryNormals(b, normals);
            for (int q = 0; q < normals[0].size(); ++q)
                valid &= std::abs(normals[0][q] - 1) < 1e-10;
            for (int i = 0; i < dim; ++i)
                valid &= m_fields[i]->GetBndConditions()[b]->GetBoundaryConditionType() == SpatialDomains::eNeumann &&
                         m_fields[i]->GetBndConditions()[b]->GetUserDefined() == tag;
            for (int i = dim; i < m_nConvectiveFields; ++i)
                valid &= m_fields[i]->GetBndConditions()[b]->GetBoundaryConditionType() == SpatialDomains::eNeumann;
        }
        else if (tag == "AWFOriginal" || tag == "AWFEquilibrium")
        {
            valid &= wallModel == tag;
            terrainBoundary = true;
            for (int i = 0; i < dim - (m_terrainWall ? 0 : 1); ++i)
                valid &= m_fields[i]->GetBndConditions()[b]->GetBoundaryConditionType() == SpatialDomains::eRobin &&
                         m_fields[i]->GetBndConditions()[b]->GetUserDefined() == tag;
            valid &= (m_terrainWall || m_fields[dim - 1]->GetBndConditions()[b]->GetBoundaryConditionType() == SpatialDomains::eDirichlet) &&
                     m_fields[m_k]->GetBndConditions()[b]->GetBoundaryConditionType() == SpatialDomains::eRobin &&
                     m_fields[m_k]->GetBndConditions()[b]->GetUserDefined() == tag &&
                     m_fields[m_epsilon]->GetBndConditions()[b]->GetBoundaryConditionType() == SpatialDomains::eDirichlet;
        }
        else if (tag == "RANSLogLaw")
        {
            valid &= m_model == RANSModel::StandardKEpsilon && wallModel == "None";
            terrainBoundary = true;
            for (int i = 0; i < dim - (m_terrainWall ? 0 : 1); ++i)
            {
                valid &= m_fields[i]->GetBndConditions()[b]->GetBoundaryConditionType() ==
                             (m_terrainWall ? SpatialDomains::eRobin : SpatialDomains::eNeumann) &&
                         m_fields[i]->GetBndConditions()[b]->GetUserDefined() == tag;
            }
            valid &= (m_terrainWall || m_fields[dim - 1]->GetBndConditions()[b]->GetBoundaryConditionType() == SpatialDomains::eDirichlet) &&
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
    if (m_terrainWall)
    {
        bool sameVelocityMaps = true;
        for (int b = 0; b < m_fields[0]->GetBndConditions().size(); ++b)
            for (int i = 1; i < dim; ++i)
                sameVelocityMaps &= m_fields[i]->GetBndConditions()[b]->GetBoundaryConditionType() ==
                                    m_fields[0]->GetBndConditions()[b]->GetBoundaryConditionType();
        RequireEverywhere(sameVelocityMaps,"terrain vector solve needs the same Cartesian velocity boundary type on each region");
    }
    int outflow = m_outflows.empty() ? 0 : 1;
    m_session->GetComm()->AllReduce(outflow, LibUtilities::ReduceMax);
    m_hasOutflow = outflow;
    m_wallModel = m_terrainWall && wallModel == "None" ? "RANSLogLaw" : wallModel;
    int hasTerrainBoundary = terrainBoundary ? 1 : 0;
    m_session->GetComm()->AllReduce(hasTerrainBoundary,LibUtilities::ReduceMax);
    RequireEverywhere(!m_terrainWall || hasTerrainBoundary,"terrain geometry requires a supported tagged wall");
    if (wallModel != "None" || m_terrainWall)
    {
        for (const auto &force : m_forcing)
            ASSERTL0(std::dynamic_pointer_cast<ForcingBody>(force),
                "The initial AWF adapter supports Body/Field forcing only, since its current forcing is evaluated within the implicit iteration");
        if (m_terrainWall)
            m_wallStrip = std::make_shared<RANSTerrainWall>(m_session,m_fields,dim,m_wallGeometry,*m_reference);
        else m_wallStrip = std::make_shared<RANSAnalyticalWallStrip>(
            m_session, m_fields, dim, m_wallOrigin, *m_reference);
    }
    // Pressure's Poisson map has an artificial pinned vertex. Mass inverses
    // used by the preconditioner must span the FULL pressure space, including
    // that vertex: a mass operator is nonsingular without any pin.
    m_pressureMass = MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(
        m_session, m_fields[0]->GetGraph(), "p", true, false);
    m_pressurePoisson = m_pressure;
    std::string linear;
    m_session->LoadSolverInfo("RANSLinearSolver", linear, "Native");
    ASSERTL0(linear == "Native" || linear == "MatrixFree", "RANSLinearSolver must be Native or MatrixFree");
    ASSERTL0(!m_terrainWall || linear == "MatrixFree","terrain wall geometry requires RANSLinearSolver=MatrixFree");
    if (linear == "MatrixFree")
    {
        double tolerance;
        m_session->LoadParameter("IterativeSolverTolerance", tolerance, 1e-13);
        // The native static-condensation inverse is a preconditioner for the
        // current-coefficient operator. Its accuracy need not be the stopping
        // tolerance of that outer solve. Keep the historical default unless
        // the independent current-operator tolerance is explicitly supplied.
        m_session->LoadParameter("RANSMatrixFreeTolerance", tolerance, tolerance);
        RANSMatrixFree::RobinOperator robin;
        if (m_wallStrip) robin = [this](int f, const Fields &v, Array<OneD, NekDouble> &r) {
            m_wallStrip->AddRobinAuditResidual(f, v, r);
        };
        bool supported = true;
        for (int f = 0; f < m_nConvectiveFields; ++f)
            for (const auto &bc : m_fields[f]->GetBndConditions())
                if (bc->GetBoundaryConditionType() == SpatialDomains::eRobin)
                    supported &= m_wallStrip && bc->GetUserDefined() == m_wallModel;
        RequireEverywhere(supported, "matrix-free Robin conditions require the analytical strip");
        m_matrixFree = std::make_shared<RANSMatrixFree>(m_fields, dim, tolerance, robin, m_gradDiv);
    }
    m_historyScheme = std::make_shared<RANSIMEX>(m_intScheme->GetOrder());
    m_intScheme = m_historyScheme;
    ASSERTL0(m_checktime == 0, "RANS constant-step history requires checkpoint intervals in steps");
}

std::string RANSVelocityCorrectionScheme::RestartSignature() const
{
    std::ostringstream s;
    s << std::setprecision(17) << "RANS-1:" << static_cast<int>(m_model) << ':'
      << m_parameters.cMu << ':' << m_parameters.cEpsilon1 << ':' << m_parameters.cEpsilon2
      << ':' << m_parameters.sigmaK << ':' << m_parameters.sigmaEpsilon << ':'
      << m_parameters.yapCoefficient << ':' << m_parameters.yapLength << ':'
      << m_kinvis << ':' << m_resolvedNu << ':' << m_wallOrigin << ':' << m_z0
      << ':' << m_heightShift << ':' << m_minimumK << ':' << m_minimumEpsilon
      << ':' << m_hasOutflow << ":pi=p/rho+2k/3";
    for (const auto &name : m_session->GetVariables()) s << ':' << name;
    if (m_wallStrip) s << ':' << m_wallStrip->RestartSignature();
    if (m_terrainWall) s << ":Newton=" << m_useNewton << ":difference=" << m_newtonDifference << ":Krylov=" << m_newtonKrylovDimension;
    if (m_preconditionedNewton) s << ":NewtonMap=PreconditionedResidual";
    if (m_session->DefinesParameter("RANSMatrixFreeTolerance"))
        s << ":currentOperatorTolerance=" << m_session->GetParameter("RANSMatrixFreeTolerance");
    if (m_startupDuration > 0) s << ":startupDuration=" << m_startupDuration;
    if (m_gradDiv > 0) s << ":gradDiv=" << m_gradDiv;
    if (m_implicitAdvection) s << ":advection=Implicit";
    s << RANSExperimentalInletSignature(m_session);
    return s.str();
}

bool RANSVelocityCorrectionScheme::CheckState(const Fields &state, bool requirePositive) const
{
    bool valid = true;
    double minK = std::numeric_limits<double>::max(), minE = minK;
    for (int f = 0; f < state.size(); ++f)
        for (int q = 0; q < state[f].size(); ++q)
            valid &= std::isfinite(state[f][q]);
    for (int q = 0; q < state[m_k].size(); ++q)
    {
        minK = std::min(minK, state[m_k][q]);
        minE = std::min(minE, state[m_epsilon][q]);
    }
    m_session->GetComm()->AllReduce(minK, LibUtilities::ReduceMin);
    m_session->GetComm()->AllReduce(minE, LibUtilities::ReduceMin);
    const bool positive = minK > m_minimumK && minE > m_minimumEpsilon;
    std::ostringstream message;
    message << std::setprecision(17)
            << (requirePositive ? "nonfinite state or k/eps positivity guard activated; min k=" :
                                  "nonfinite linear preconditioner image; min k=")
            << minK << ", min eps=" << minE << "; no floor is applied";
    RequireEverywhere(valid && (!requirePositive || positive), message.str());
    return positive;
}

void RANSVelocityCorrectionScheme::v_SetInitialConditions(NekDouble time, bool dump, int domain)
{
    EquationSystem::v_SetInitialConditions(time, dump, domain);
    // Bind history to the imported field before VCS reapplies its startup
    // Dirichlet values. The dynamic AWF EH value is restored from history.
    m_importedCoefficients = Fields(m_fields.size());
    for (int f = 0; f < m_fields.size(); ++f)
    {
        m_importedCoefficients[f] = Array<OneD, NekDouble>(m_fields[f]->GetNcoeffs());
        Vmath::Vcopy(m_fields[f]->GetNcoeffs(),m_fields[f]->GetCoeffs(),1,m_importedCoefficients[f],1);
    }
}

void RANSVelocityCorrectionScheme::v_DoInitialise(bool dump)
{
    std::string mode, input;
    m_session->LoadSolverInfo("RANSRestartMode", mode, "Auto");
    ASSERTL0(mode == "Auto" || mode == "History" || mode == "Cold",
             "RANSRestartMode must be Auto, History or Cold");
    int fileVariables = 0;
    for (const auto &name : m_session->GetVariables())
        if (m_session->DefinesFunction("InitialConditions", name) &&
            m_session->GetFunctionType("InitialConditions", name) == LibUtilities::eFunctionTypeFile)
        {
            const auto file = m_session->GetFunctionFilename("InitialConditions", name);
            ASSERTL0(input.empty() || input == file, "RANS restart fields must use one common file");
            input = file; ++fileVariables;
        }
    LibUtilities::FieldMetaDataMap stored;
    if (fileVariables)
    {
        auto path = std::filesystem::path(input);
        if (std::filesystem::is_directory(path)) path /= "Info.xml";
        auto reader = LibUtilities::FieldIO::CreateForFile(m_session, path.string());
        reader->ImportFieldMetaData(path.string(), stored);
        if (stored.count("RANSConfiguration"))
            ASSERTL0(stored.at("RANSConfiguration") == RestartSignature(),
                     "RANS restart model/configuration mismatch");
    }
    const bool restore = mode != "Cold" && stored.count("RANSHistoryFile");
    ASSERTL0(mode != "History" || restore, "RANS history restart requires history metadata");
    ASSERTL0(!restore || (fileVariables == m_fields.size() && stored.count("Time") &&
             stored.count("RANSConfiguration")), "RANS history restart requires all primary fields and metadata");
    VCSWeakPressure::v_DoInitialise(dump);
    ASSERTL0(m_startupDuration == 0 || m_time == 0 || m_time >= m_startupDuration,
        "startup weak load is not checkpointed; restart from a post-startup checkpoint");
    if (restore)
    {
        std::string path;
        m_session->LoadSolverInfo("RANSHistoryPath", path, stored.at("RANSHistoryFile"));
        m_historyStepOffset = RANSHistory::Read(path, m_fields, m_importedCoefficients, *m_historyScheme,
            RestartSignature(), m_time, m_timestep);
    }
    else if (fileVariables && m_session->GetComm()->GetRank() == 0)
        std::cout << "RANS cold restart: native IMEX startup; no saved multistep history restored\n";
    m_fieldMetaDataMap.erase("RANSHistoryFile");
    m_importedCoefficients = Fields();
    m_fieldMetaDataMap["RANSConfiguration"] = RestartSignature();
    Fields initial(m_nConvectiveFields);
    for (int i = 0; i < initial.size(); ++i)
        initial[i] = m_fields[i]->UpdatePhys();
    CheckState(initial);
    if (m_wallStrip) m_wallStrip->CheckInitialState(initial);
    int auditStage = 0;
    m_session->LoadParameter("RANSStageAudit", auditStage, 0);
    ASSERTL0(auditStage == 0 || auditStage == 1, "RANSStageAudit must be 0 or 1");
    m_auditStage = auditStage;
    int monitor = 0;
    m_session->LoadParameter("RANSReferenceMonitor", monitor, 0);
    ASSERTL0(monitor == 0 || monitor == 1, "RANSReferenceMonitor must be 0 or 1");
    m_monitorReference = monitor;
    m_fieldMetaDataMap["RANSModel"] = m_model == RANSModel::CLSHighRe ? "CLSHighRe" : "StandardKEpsilon";
    m_fieldMetaDataMap["RANSDissipation"] = "high-Re: eps_total=eps_transport; molecular mode adds diffusion only";
    m_fieldMetaDataMap["RANSPhysicalFrame"] = m_velocity.size() == 2 ? "mesh(x,y)=physical(x,z)" : "physical(x,y,z)";
    m_fieldMetaDataMap["RANSWallModel"] = m_wallStrip ? m_wallModel + " conservative strip" : "None or standard log law";
    m_fieldMetaDataMap["RANSPressure"] = "pi=p/rho+2k/3";
    m_fieldMetaDataMap["RANSResolvedNu"] = std::to_string(m_resolvedNu);
    if (m_gradDiv > 0) m_fieldMetaDataMap["RANSGradDiv"] = std::to_string(m_gradDiv);
    if (m_startupDuration > 0)
        m_fieldMetaDataMap["RANSStartupDuration"] = std::to_string(m_startupDuration);
    if (m_session->GetComm()->GetRank() == 0)
    {
        if (m_monitorReference)
        {
            m_referenceDiagnostics.open(m_session->GetSessionName() + ".rans-errors.csv");
            ASSERTL0(m_referenceDiagnostics.good(), "Cannot open RANS reference monitor");
            m_referenceDiagnostics << std::setprecision(17)
                << "time,field,l2_error,reference_l2,normalized_l2,pointwise_relative_linf,first_layer_relative_l2,upper_layer_relative_l2\n";
        }
        if (m_auditStage)
        {
            m_residualDiagnostics.open(m_session->GetSessionName() + ".rans-residual.csv");
            ASSERTL0(m_residualDiagnostics.good(), "Cannot open RANS residual audit");
            m_residualDiagnostics << "time,phase,field,free_coefficient_l2,free_coefficient_linf\n" << std::setprecision(17);
        }
        m_diagnostics.open(m_session->GetSessionName() + ".rans.csv");
        ASSERTL0(m_diagnostics.good(), "Cannot open RANS diagnostics");
        m_diagnostics << "time,iterations,change,pressure_correction,divergence_l2,weak_divergence_l2,min_k,min_eps,min_nut,max_nut,live_global_systems,guard_activations,matrix_free_iterations,newton_linear_iterations,newton_backtracks,startup_weight,nonpositive_linear_images\n" << std::setprecision(17);
    }
    if (m_monitorReference) MonitorReference();
}

void RANSVelocityCorrectionScheme::v_GenerateSummary(SummaryList &summary)
{
    VCSWeakPressure::v_GenerateSummary(summary);
    AddSummaryItem(summary, "RANS model", m_model == RANSModel::CLSHighRe ?
        "High-Re CLS, full cubic stress and algebraic Yap" : "Standard k-epsilon, full symmetric stress");
    AddSummaryItem(summary, "RANS resolved molecular viscosity", m_resolvedNu);
    AddSummaryItem(summary, "RANS neutral kappa", m_reference->Kappa());
    AddSummaryItem(summary, "RANS wall model", m_wallModel);
    AddSummaryItem(summary, "RANS implicit stage", m_useNewton ? "Newton-Krylov coupled stage" : "nonlinear Picard + pressure defect iteration");
    if (m_useNewton) AddSummaryItem(summary,"RANS Newton map",m_preconditionedNewton ?
        "current residual with frozen-inverse sweeps" : "current Helmholtz solves");
    AddSummaryItem(summary, "RANS pressure correction tolerance", m_pressureTolerance);
    AddSummaryItem(summary, "RANS advection treatment", m_implicitAdvection ? "implicit nonlinear stage" : "explicit IMEX history");
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
    m_fieldMetaDataMap["RANSObservationConfiguration"] = RANSObservationConfiguration(m_session);
    m_fieldMetaDataMap.erase("RANSInitializationOnly");
    const int dim = m_velocity.size(), nq = m_fields[0]->GetTotPoints();
    Gradients gradient(dim);
    Fields state(m_nConvectiveFields);
    for (int i = 0; i < state.size(); ++i)
    {
        state[i] = Array<OneD, NekDouble>(nq);
        m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), state[i]);
    }
    CheckState(state);
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
    std::string outputError;
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
    catch (const std::exception &error) { valid = false; outputError = error.what(); }
    RequireEverywhere(valid, "cannot output an inadmissible RANS state: "+outputError);
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
    if (m_implicitAdvection)
        for (auto &field : output) Vmath::Zero(field.size(),field,1);
    else EvaluateAdvectionTerms(input, output, time);
    for (auto &force : m_forcing)
        force->Apply(m_fields, input, output, time);
    // No legacy curl-curl pressure extrapolation. The implicit stage enforces
    // continuity with the actual variable symmetric stress and stage BCs.
}

void RANSVelocityCorrectionScheme::UpdateLogWall(const Fields &state)
{
    if (m_terrainWall) return; // vector terrain adapter supplies this wall policy
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
    Fields &gain, Fields &loss, NekDouble time)
{
    const int nq = m_fields[0]->GetTotPoints(), dim = m_velocity.size();
    CheckState(state);
    for (int i = 0; i < dim; ++i)
        for (int j = 0; j < dim; ++j)
            m_fields[i]->PhysDeriv(DirCartesianMap[j], state[i], gradient[i][j]);
    if (m_cmuBranches.size() != nq) m_cmuBranches.resize(nq);
    bool valid = true;
    try
    {
      for (int q = 0; q < nq; ++q)
      {
        RANSTensor g{};
        const auto axis = [dim](int i) { return dim == 2 && i == 1 ? 2 : i; };
        for (int i = 0; i < dim; ++i)
            for (int j = 0; j < dim; ++j) g[axis(i)][axis(j)] = gradient[i][j][q];
        const auto v = m_freezeCmuBranch ? EvaluateRANSClosureLinearization(m_model,
            state[m_k][q],state[m_epsilon][q],g,m_wallDistance[q],m_parameters,m_resolvedNu,m_cmuBranches[q]) :
            EvaluateRANSClosure(m_model,state[m_k][q],state[m_epsilon][q],g,m_wallDistance[q],m_parameters,m_resolvedNu);
        if (!m_freezeCmuBranch) m_cmuBranches[q] = v.strainInvariant >= v.rotationInvariant;
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
    if (m_implicitAdvection)
    {
        Fields advection(state.size());
        for (auto &field : advection) field = Array<OneD,NekDouble>(nq,0.0);
        EvaluateAdvectionTerms(state,advection,time);
        for (int i = 0; i < state.size(); ++i)
        {
            // The native advection result is -U.grad(phi). Its entire
            // derivative is included by the coupled nonlinear stage.
            if (i < dim) Vmath::Vcopy(nq,advection[i],1,gain[i],1);
            else Vmath::Vadd(nq,advection[i],1,gain[i],1,gain[i],1);
        }
    }
}

void RANSVelocityCorrectionScheme::UpdateWall(
    const Fields &state, const Fields &predictor, const Gradients &gradient,
    const Array<OneD, NekDouble> &pressure, NekDouble time, NekDouble alpha)
{
    if (!m_wallStrip) return;
    const int nf = state.size(), dim = m_velocity.size(), nq = pressure.size();
    Fields advection(nf), body(nf), pressureForce(dim), convection(dim);
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
        convection[i] = Array<OneD, NekDouble>(nq);
        m_fields[0]->PhysDeriv(DirCartesianMap[i], mechanical, pressureForce[i]);
        Vmath::Smul(nq, -1.0, advection[i], 1, convection[i], 1);
    }
    if (m_implicitAdvection)
    {
        // The strip's inertial predictor normally includes explicit
        // convection. Supply its CURRENT contribution instead, including k,
        // so moving bulk transport does not omit the strip inventory flux.
        Fields adjusted(nf);
        for (int i = 0; i < nf; ++i)
        {
            adjusted[i] = Array<OneD,NekDouble>(nq);
            Vmath::Svtvp(nq,1/alpha,advection[i],1,predictor[i],1,adjusted[i],1);
        }
        m_wallStrip->Update(state, adjusted, gradient, pressureForce, convection, body, alpha);
    }
    else m_wallStrip->Update(state, predictor, gradient, pressureForce, convection, body, alpha);
}

void RANSVelocityCorrectionScheme::UpdateOutflow(
    const Array<OneD, NekDouble> &pressure, NekDouble time)
{
    const int dim = m_velocity.size();
    for (int b : m_outflows)
    {
        Array<OneD, NekDouble> trace;
        m_pressure->ExtractPhysToBnd(b, pressure, trace);
        auto boundary = m_fields[0]->GetBndCondExpansions()[b];
        Array<OneD, NekDouble> x(trace.size()), y(trace.size()), z(trace.size()), values(trace.size());
        boundary->GetCoords(x, y, z);
        Fields normal;
        m_fields[0]->GetBoundaryNormals(b, normal);
        // XML prescribes total traction (T - pi I)n. Native Helmholtz
        // expects Tn, so the CURRENT pi trace completes its Neumann load.
        // All pressure test functions remain free: continuity determines
        // the outlet flux and the traction determines the pressure gauge.
        for (int i = 0; i < dim; ++i)
        {
            const auto bc = m_fields[i]->GetBndConditions()[b];
            bc->GetEquation()->Evaluate(x, y, z, time, values);
            Vmath::Vvtvp(trace.size(), trace, 1, normal[i], 1, values, 1, values, 1);
            boundary = m_fields[i]->GetBndCondExpansions()[b];
            boundary->IProductWRTBase(values, boundary->UpdateCoeffs());
        }
    }
}

void RANSVelocityCorrectionScheme::MonitorReference()
{
    const int nq = m_fields[0]->GetTotPoints();
    std::vector<std::string> names = m_session->GetVariables();
    Fields values(names.size());
    for (int f = 0; f < values.size(); ++f)
    {
        values[f] = Array<OneD, NekDouble>(nq);
        m_fields[f]->BwdTrans(m_fields[f]->GetCoeffs(), values[f]);
    }
    std::vector<Array<OneD, NekDouble>> extra;
    std::vector<std::string> extraNames;
    v_ExtraFldOutput(extra, extraNames);
    const int primary = values.size();
    Fields all(primary + extra.size());
    for (int i = 0; i < primary; ++i) all[i] = values[i];
    for (int i = 0; i < extra.size(); ++i)
    {
        all[primary+i] = Array<OneD, NekDouble>(nq);
        m_fields[0]->BwdTrans(extra[i], all[primary+i]);
        names.push_back(extraNames[i]);
    }
    Array<OneD, NekDouble> x(nq, 0.0), y(nq, 0.0), z(nq, 0.0), reference(nq), error(nq), work(nq), ones(nq, 1.0);
    m_fields[0]->GetCoords(x, y, z);
    double star = 1, nextHeight = 4;
    m_session->LoadParameter("Ustar", star, 1.0);
    m_session->LoadParameter("RANSStripPointN", nextHeight, 4.0);
    const double volume = m_fields[0]->Integral(ones);
    std::array<double, 2> layerVolume;
    for (int layer = 0; layer < 2; ++layer)
    {
        for (int q = 0; q < nq; ++q) work[q] = ((m_wallDistance[q] <= nextHeight) == (layer == 0)) ? 1.0 : 0.0;
        layerVolume[layer] = m_fields[0]->Integral(work);
    }
    ASSERTL0(std::isfinite(star) && star > 0, "Reference monitor requires positive Ustar");
    for (int f = 0; f < all.size(); ++f)
    {
        if (!m_session->DefinesFunction("ExactSolution", names[f])) continue;
        m_session->GetFunction("ExactSolution", names[f])->Evaluate(x, y, z, m_time, reference);
        Vmath::Vsub(nq, all[f], 1, reference, 1, error, 1);
        const double norm = m_fields[0]->L2(reference), err = m_fields[0]->L2(error);
        const double fallback = f < m_velocity.size() ? star : star*star;
        double relativeMax = 0;
        for (int q = 0; q < nq; ++q)
            relativeMax = std::max(relativeMax, std::abs(error[q]) /
                (std::abs(reference[q]) > 1e-30 ? std::abs(reference[q]) : fallback));
        relativeMax = Maximum(relativeMax);
        std::array<double, 2> layerError, layerReference;
        for (int layer = 0; layer < 2; ++layer)
        {
            for (int q = 0; q < nq; ++q)
                work[q] = ((m_wallDistance[q] <= nextHeight) == (layer == 0)) ? error[q]*error[q] : 0.0;
            layerError[layer] = std::sqrt(std::max(0.0, m_fields[0]->Integral(work)));
            for (int q = 0; q < nq; ++q)
                work[q] = ((m_wallDistance[q] <= nextHeight) == (layer == 0)) ? reference[q]*reference[q] : 0.0;
            layerReference[layer] = std::sqrt(std::max(0.0, m_fields[0]->Integral(work)));
        }
        if (m_referenceDiagnostics.is_open())
            m_referenceDiagnostics << m_time << ',' << names[f] << ',' << err << ',' << norm << ','
                << err / (norm > 1e-30 ? norm : fallback*std::sqrt(volume)) << ',' << relativeMax << ','
                << layerError[0]/(layerReference[0] > 1e-30 ? layerReference[0] : fallback*std::sqrt(std::max(layerVolume[0], 1e-30))) << ','
                << layerError[1]/(layerReference[1] > 1e-30 ? layerReference[1] : fallback*std::sqrt(std::max(layerVolume[1], 1e-30))) << '\n' << std::flush;
    }
}

bool RANSVelocityCorrectionScheme::v_PostIntegrate(int step)
{
    const bool stop = VelocityCorrectionScheme::v_PostIntegrate(step);
    if (m_monitorReference) MonitorReference();
    if ((m_checksteps && !((step+1)%m_checksteps)) || stop ||
        step+1 == m_steps || (m_fintime > 0 && m_time >= m_fintime-1e-12))
    {
        const auto completed = m_historyStepOffset + step + 1;
        const auto directory = std::filesystem::absolute(m_session->GetSessionName() +
            ".rans-" + std::to_string(completed)).string();
        const auto published = RANSHistory::Write(directory, m_fields, *m_historyScheme,
            RestartSignature(), completed, m_time, m_timestep);
        m_fieldMetaDataMap["RANSHistoryFile"] = published;
    }
    return stop;
}

void RANSVelocityCorrectionScheme::AuditStage(
    const Fields &state, const Fields &predictor, const Fields &diffusion,
    const Fields &gain, const Fields &loss,
    const Array<OneD, NekDouble> &pressure, NekDouble time, NekDouble alpha, int phase)
{
    // Assemble the complete nonlinear stage residual on the actual trial
    // state. No solve, equilibrium subtraction or analytical derivative enters
    // this audit. Norms are in the native free global coefficient basis.
    const int dim = m_velocity.size(), nq = pressure.size();
    Array<OneD,NekDouble> divergence(nq,0.0),derivative(nq);
    if (m_gradDiv > 0)
        for (int i = 0; i < dim; ++i)
        {
            m_fields[i]->PhysDeriv(DirCartesianMap[i],state[i],derivative);
            Vmath::Vadd(nq,derivative,1,divergence,1,divergence,1);
        }
    auto record = [&](const MultiRegions::ExpListSharedPtr &field,
                      const Array<OneD, NekDouble> &local, const std::string &name)
    {
        auto cg = std::dynamic_pointer_cast<MultiRegions::ContField>(field);
        const auto map = cg->GetLocalToGlobalMap();
        Array<OneD, NekDouble> global(map->GetNumGlobalCoeffs());
        cg->Assemble(local, global);
        double squared = 0, maximum = 0;
        for (int j = map->GetNumGlobalDirBndCoeffs(); j < global.size(); ++j)
            if (map->GetGlobalToUniversalMapUnique(j) > 0)
            {
                squared += global[j] * global[j];
                maximum = std::max(maximum, std::abs(global[j]));
            }
        m_session->GetComm()->AllReduce(squared, LibUtilities::ReduceSum);
        maximum = Maximum(maximum);
        if (m_residualDiagnostics.is_open())
            m_residualDiagnostics << time << ',' << phase << ',' << name << ','
                << std::sqrt(squared) << ',' << maximum << '\n' << std::flush;
    };
    for (int i = 0; i < state.size(); ++i)
    {
        bool supported = true;
        for (const auto &bc : m_fields[i]->GetBndConditions())
            if (bc->GetBoundaryConditionType() == SpatialDomains::eRobin)
                supported &= m_wallStrip && bc->GetUserDefined() == m_wallModel;
        RequireEverywhere(supported, "RANSStageAudit supports only the analytical strip's Robin mass");
        const int nc = m_fields[i]->GetNcoeffs();
        Array<OneD, NekDouble> source(nq), local(nc), weak(nc), boundary(nc, 0.0);
        for (int q = 0; q < nq; ++q)
            source[q] = alpha * (state[i][q] - predictor[i][q]) - gain[i][q] + loss[i][q] * state[i][q];
        if (i < dim)
        {
            Array<OneD, NekDouble> gradient(nq);
            m_pressure->PhysDeriv(DirCartesianMap[i], pressure, gradient);
            Vmath::Vadd(nq, source, 1, gradient, 1, source, 1);
        }
        m_fields[i]->IProductWRTBase(source, local);
        Fields flux(dim);
        for (int j = 0; j < dim; ++j)
        {
            flux[j] = Array<OneD, NekDouble>(nq);
            m_fields[i]->PhysDeriv(DirCartesianMap[j], state[i], flux[j]);
            if (i < dim)
            {
                Array<OneD, NekDouble> transpose(nq);
                m_fields[j]->PhysDeriv(DirCartesianMap[i], state[j], transpose);
                Vmath::Vadd(nq, flux[j], 1, transpose, 1, flux[j], 1);
            }
            Vmath::Vmul(nq, diffusion[i], 1, flux[j], 1, flux[j], 1);
            if (i < dim) Vmath::Vsub(nq, flux[j], 1, m_nonlinearStress[i][j], 1, flux[j], 1);
            if (m_gradDiv > 0 && i < dim && i == j)
                Vmath::Svtvp(nq,m_gradDiv,divergence,1,flux[j],1,flux[j],1);
        }
        m_fields[i]->IProductWRTDerivBase(flux, weak);
        Vmath::Vadd(nc, local, 1, weak, 1, local, 1);
        m_fields[i]->ImposeNeumannConditions(boundary);
        m_fields[i]->ImposeRobinConditions(boundary);
        Vmath::Vsub(nc, local, 1, boundary, 1, local, 1);
        if (m_wallStrip) m_wallStrip->AddRobinAuditResidual(i, state, local);
        if (phase == -1)
            m_startupLoad[i] = local;
        else
        {
            if (m_startupLoad.size() && StartupWeight(time) > 0)
                Vmath::Svtvp(nc, -StartupWeight(time), m_startupLoad[i], 1,
                             local, 1, local, 1);
            record(m_fields[i], local, m_session->GetVariable(i));
        }
    }
    if (phase == -1) return;
    Fields velocity(dim);
    for (int i = 0; i < dim; ++i) velocity[i] = state[i];
    Array<OneD, NekDouble> weak(m_pressureMass->GetNcoeffs());
    m_pressureMass->IProductWRTDerivBase(velocity, weak);
    Vmath::Smul(weak.size(), -1.0, weak, 1, weak, 1);
    auto cg = std::dynamic_pointer_cast<MultiRegions::ContField>(m_pressureMass);
    const auto map = cg->GetLocalToGlobalMap();
    int offset = 0;
    for (int b = 0; b < m_pressureMass->GetBndCondExpansions().size(); ++b)
    {
        auto boundary = m_pressureMass->GetBndCondExpansions()[b];
        if (m_pressureMass->GetBndConditions()[b]->GetBoundaryConditionType() != SpatialDomains::ePeriodic)
        {
            Fields trace(dim);
            for (int i = 0; i < dim; ++i) m_fields[i]->ExtractPhysToBnd(b, state[i], trace[i]);
            Array<OneD, NekDouble> coefficients(boundary->GetNcoeffs());
            boundary->NormVectorIProductWRTBase(trace, coefficients);
            for (int j = 0; j < coefficients.size(); ++j)
                weak[map->GetBndCondCoeffsToLocalCoeffsMap()[offset+j]] +=
                    (map->GetSignChange() ? map->GetBndCondCoeffsToLocalCoeffsSign()[offset+j] : 1.0) * coefficients[j];
        }
        offset += boundary->GetNcoeffs();
    }
    record(m_pressureMass, weak, "continuity");
    if (m_wallStrip && phase == 0) m_wallStrip->WritePointDiagnostics(time, phase);
}

void RANSVelocityCorrectionScheme::SolveField(
    int field, const Array<OneD, NekDouble> &rhs,
    const Array<OneD, NekDouble> &diffusion,
    const Array<OneD, NekDouble> &mass, Array<OneD, NekDouble> &output, NekDouble alpha,
    const Array<OneD, NekDouble> *newtonTrial)
{
    if (m_matrixFree)
    {
        m_linearIterations += m_matrixFree->Solve(field, rhs, diffusion, mass, alpha, output,newtonTrial);
        return;
    }
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
    if (m_hasOutflow && m_pressureInverseAlpha != alpha)
    {
        Array<OneD, NekDouble> one(nq,1.0);
        const NekDouble meanMu = m_pressureMass->Integral(mu)/m_pressureMass->Integral(one);
        ASSERTL0(std::isfinite(meanMu) && meanMu > 0,"Invalid pressure inverse viscosity");
        // Resolving the velocity boundary layer gives an outlet pressure
        // impedance O(sqrt(alpha/mu)). Freeze it only in this preconditioner
        // for each IMEX inertia value; all physical coefficients stay current.
        m_pressurePoisson = std::make_shared<RANSOutletPressureInverse>(
            m_session,m_fields[0]->GetGraph(),m_outflows,std::sqrt(alpha/meanMu));
        m_pressureInverseAlpha = alpha;
    }
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
    // Assemble before normalizing: local volume/boundary terms can be
    // large even when their assembled divergence is tiny. Redistribute the
    // complete RHS so native solver tolerance has no earlier-RHS dependence.
    auto normalizedLoad = [&](const MultiRegions::ExpListSharedPtr &field,
                              const Array<OneD,NekDouble> &load,Array<OneD,NekDouble> &distributed) {
        const auto map = std::dynamic_pointer_cast<MultiRegions::ContField>(field)->GetLocalToGlobalMap();
        Array<OneD,NekDouble> global(map->GetNumGlobalCoeffs()),counts(global.size()),local(nc);
        map->Assemble(load,global);
        for (int i = 0; i < nc; ++i) local[i] = map->GetLocalToGlobalSign(i);
        map->Assemble(local,counts);
        double norm = 0;
        for (int i = map->GetNumGlobalDirBndCoeffs(); i < global.size(); ++i)
            if (map->GetGlobalToUniversalMapUnique(i) > 0) norm += global[i]*global[i];
        m_session->GetComm()->AllReduce(norm,LibUtilities::ReduceSum);norm = std::sqrt(norm);
        map->GlobalToLocal(global,distributed);
        for (int i = 0; i < nc; ++i) distributed[i] /= counts[map->GetLocalToGlobalMap(i)]*(norm > 0 ? norm : 1);
        return norm;
    };
    if (m_useNewton)
    {
        const double norm = normalizedLoad(m_pressurePoisson,weak,rhs);
        std::vector<std::pair<MultiRegions::ExpListSharedPtr,Array<OneD,NekDouble>>> saved;
        for (int b = 0; b < m_pressurePoisson->GetBndCondExpansions().size(); ++b)
        {
            if (m_pressurePoisson->GetBndConditions()[b]->GetBoundaryConditionType() == SpatialDomains::ePeriodic) continue;
            auto boundary = m_pressurePoisson->GetBndCondExpansions()[b];
            Array<OneD,NekDouble> copy(boundary->GetNcoeffs());
            Vmath::Vcopy(copy.size(),boundary->GetCoeffs(),1,copy,1);saved.emplace_back(boundary,copy);
            Vmath::Zero(copy.size(),boundary->UpdateCoeffs(),1);
        }
        Vmath::Zero(nc,m_pressurePoisson->UpdateCoeffs(),1);
        Vmath::Zero(nq,increment,1);
        if (norm > 0)
        {
            StdRegions::ConstFactorMap factors;factors[StdRegions::eFactorLambda] = 0;
            m_pressurePoisson->HelmSolve(rhs,m_pressurePoisson->UpdateCoeffs(),factors,
                StdRegions::NullVarCoeffMap,StdRegions::NullVarFactorsMap,NullNekDouble1DArray,false);
            m_pressurePoisson->BwdTrans(m_pressurePoisson->GetCoeffs(),increment);
            Vmath::Smul(nq,alpha*norm,increment,1,increment,1);
        }
        for (const auto &entry : saved)
            Vmath::Vcopy(entry.second.size(),entry.second,1,entry.first->UpdateCoeffs(),1);
    }
    else
    {
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
        if (m_hasOutflow)
            for (int b = 0; b < m_pressurePoisson->GetBndCondExpansions().size(); ++b)
            {
                if (m_pressurePoisson->GetBndConditions()[b]->GetBoundaryConditionType() == SpatialDomains::ePeriodic) continue;
                auto to = m_pressurePoisson->GetBndCondExpansions()[b];
                Vmath::Vcopy(to->GetNcoeffs(),m_pressure->GetBndCondExpansions()[b]->GetCoeffs(),1,to->UpdateCoeffs(),1);
            }
        // Homogeneous initial guess. Neither preconditioner prescribes the final
        // pressure multiplier at an open outlet.
        Vmath::Zero(nc, m_pressurePoisson->UpdateCoeffs(), 1);
        m_pressurePoisson->HelmSolve(rhs, m_pressurePoisson->UpdateCoeffs(), factors,
            StdRegions::NullVarCoeffMap, StdRegions::NullVarFactorsMap,
            NullNekDouble1DArray, false);
        m_pressurePoisson->BwdTrans(m_pressurePoisson->GetCoeffs(), increment);
    }
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
    if (m_useNewton)
    {
        const double norm = normalizedLoad(m_pressureMass,weak,rhs);
        Vmath::Zero(nc,projection,1);
        if (norm > 0) m_pressureMass->MultiplyByInvMassMatrix(rhs,projection);
        m_pressureMass->BwdTrans(projection,tmp);
        Vmath::Smul(nq,norm,tmp,1,tmp,1);
    }
    else
    {
        m_pressureMass->MultiplyByInvMassMatrix(weak, projection);
        m_pressureMass->BwdTrans(projection, tmp);
    }
    weakDivergence = m_pressureMass->L2(tmp);
    if (!m_nonlinearTrial && m_session->DefinesParameter("RANSDebug"))
    {
        const double compatibility = m_pressureMass->Integral(tmp);
        if (m_session->GetComm()->GetRank() == 0)
            std::cout << "RANS pressure net flux " << compatibility << std::endl;
    }
    for (int q = 0; q < nq; ++q) tmp[q] *= mu[q]+m_gradDiv;
    double projectionScale = 1;
    if (m_useNewton)
    {
        projectionScale = 0;
        for (double value : tmp) projectionScale = std::max(projectionScale,std::abs(value));
        projectionScale = Maximum(projectionScale);
        if (projectionScale > 0) Vmath::Smul(nq,1/projectionScale,tmp,1,tmp,1);
    }
    if (projectionScale > 0)
    {
        m_pressureMass->FwdTrans(tmp, projection);
        m_pressureMass->BwdTrans(projection, tmp);
        if (m_useNewton) Vmath::Smul(nq,projectionScale,tmp,1,tmp,1);
    }
    Vmath::Vsub(nq, increment, 1, tmp, 1, increment, 1);
    Array<OneD, NekDouble> ones(nq, 1.0);
    const double mean = m_hasOutflow ? 0.0 : m_pressure->Integral(increment) / m_pressure->Integral(ones);
    NekDouble change = 0;
    for (int q = 0; q < nq; ++q)
    {
        increment[q] -= mean;
        pressure[q] += m_pressureRelaxation * increment[q];
        change = std::max(change, std::abs(increment[q]));
    }
    return Maximum(change);
}

NekDouble RANSVelocityCorrectionScheme::StartupWeight(NekDouble time) const
{
    if (m_startupDuration == 0 || time >= m_startupDuration) return 0;
    const double r = std::max(0.0, time / m_startupDuration);
    const double s = 1-r;
    return s*s*s*(1+3*r+6*r*r);
}

void RANSVelocityCorrectionScheme::v_SolveUnsteadyStokesSystem(
    const Array<OneD, const Array<OneD, NekDouble>> &input,
    Fields &output, NekDouble time, NekDouble lambda)
{
    ASSERTL0(lambda > 0, "RANS stage requires positive implicit timestep");
    const int nf = m_nConvectiveFields, nq = m_fields[0]->GetTotPoints();
    const int nc = m_fields[0]->GetNcoeffs(), dim = m_velocity.size();
    const double alpha = 1.0 / lambda;
    m_linearIterations = 0;
    Fields state(nf), next(nf), diffusion(nf), gain(nf), loss(nf), predictor(nf);
    Gradients gradient(dim);
    for (int i = 0; i < nf; ++i)
    {
        state[i] = Array<OneD, NekDouble>(nq);
        next[i] = Array<OneD, NekDouble>(nq);
        diffusion[i] = Array<OneD, NekDouble>(nq);
        gain[i] = Array<OneD, NekDouble>(nq, 0.0);
        loss[i] = Array<OneD, NekDouble>(nq, 0.0);
        predictor[i] = input[i];
        Vmath::Vcopy(nq, m_fields[i]->GetPhys(), 1, state[i], 1);
    }
    for (int i = 0; i < dim; ++i)
    {
        gradient[i] = Fields(dim);
        for (int j = 0; j < dim; ++j) gradient[i][j] = Array<OneD, NekDouble>(nq);
    }
    Array<OneD, NekDouble> pressure(nq), rhs(nc), cross(nc), tmp(nq), mass(nq), gradP(nq);
    Fields velocityRhs(dim);
    for (auto &v : velocityRhs) v = Array<OneD,NekDouble>(nc);
    m_pressure->BwdTrans(m_pressure->GetCoeffs(), pressure);
    Array<OneD, NekDouble> ones(nq, 1.0);
    const double pressureMean = m_hasOutflow ? 0.0 : m_pressure->Integral(pressure) / m_pressure->Integral(ones);
    for (int q = 0; q < nq; ++q) pressure[q] -= pressureMean;
    SetBoundaryConditions(time);
    Array<OneD,NekDouble> volumeMetric(nq),wallMetric;
    std::vector<double> iterationMetric,iterationScales(nf+2,1.0),differenceScales(nf+2,1.0);
    if (m_useNewton)
    {
        for (int e = 0; e < m_fields[0]->GetExpSize(); ++e)
        {
            auto expansion = m_fields[0]->GetExp(e);
            Array<OneD,NekDouble> one(expansion->GetTotPoints(),1.0),out = volumeMetric+m_fields[0]->GetPhys_Offset(e);
            expansion->MultiplyByQuadratureMetric(one,out);
        }
        wallMetric = m_wallStrip->IterationMetric();
        iterationMetric.resize((nf+1)*nq+wallMetric.size());
        double maximumK = 0;
        for (double value : state[m_k]) maximumK = std::max(maximumK,value);
        maximumK = Maximum(maximumK);
        double velocityScale = std::sqrt(maximumK);
        for (int q = 0; q < nq; ++q)
        {
            double speed = 0;
            for (int i = 0; i < dim; ++i) speed = std::hypot(speed,state[i][q]);
            velocityScale = std::max(velocityScale,speed);
        }
        velocityScale = Maximum(velocityScale);
        for (int i = 0; i < nf; ++i)
        {
            double scale = 0;
            for (double value : state[i]) scale = std::max(scale,std::abs(value));
            scale = Maximum(scale);
            if (i < dim) scale = velocityScale;
            iterationScales[i] = differenceScales[i] = scale;
            for (int q = 0; q < nq; ++q) iterationMetric[i*nq+q] = volumeMetric[q]/(scale*scale);
        }
        const double volume = m_fields[0]->Integral(ones);
        double wallVolume = 0,height;
        for (double value : wallMetric) wallVolume += value;
        m_session->GetComm()->AllReduce(wallVolume,LibUtilities::ReduceSum);
        m_session->LoadParameter("RANSStripHeight",height,2);
        const double length = volume*height/wallVolume;
        Closure(state,gradient,diffusion,gain,loss,time);
        const double viscosity = m_fields[0]->Integral(diffusion[0])/volume;
        // The implicit velocity/pressure blocks balance at the unsteady
        // inertia scale, rather than the mean-flow dynamic pressure U^2.
        const double pressureScale = velocityScale*(alpha*length+viscosity/length)*
            std::sqrt(m_pressureRelaxation/m_relaxation);
        iterationScales[nf] = pressureScale;
        differenceScales[nf] = velocityScale*velocityScale;
        for (double value : pressure) differenceScales[nf] = std::max(differenceScales[nf],std::abs(value));
        differenceScales[nf] = Maximum(differenceScales[nf]);
        iterationScales[nf+1] = differenceScales[nf+1] = maximumK;
        for (int q = 0; q < nq; ++q) iterationMetric[nf*nq+q] = volumeMetric[q]/(pressureScale*pressureScale);
        for (int q = 0; q < wallMetric.size(); ++q) iterationMetric[(nf+1)*nq+q] = wallMetric[q]/(maximumK*maximumK);

    }
    auto pack = [&]() {
        const auto wall = m_wallStrip->IterationState();
        std::vector<double> result((nf+1)*nq+wall.size());
        for (int i = 0; i < nf; ++i) for (int q = 0; q < nq; ++q) result[i*nq+q] = state[i][q];
        for (int q = 0; q < nq; ++q) result[nf*nq+q] = pressure[q];
        for (int q = 0; q < wall.size(); ++q) result[(nf+1)*nq+q] = wall[q];
        return result;
    };
    auto dotIteration = [&](const std::vector<double> &a,const std::vector<double> &b) {
        double sum = 0;
        for (int i = 0; i < a.size(); ++i) sum += iterationMetric[i]*a[i]*b[i];
        m_session->GetComm()->AllReduce(sum,LibUtilities::ReduceSum);
        return sum;
    };
    auto unpack = [&](const std::vector<double> &values) {
        for (int i = 0; i < nf; ++i) for (int q = 0; q < nq; ++q) state[i][q] = values[i*nq+q];
        for (int q = 0; q < nq; ++q) pressure[q] = values[nf*nq+q];
        // The pressure iteration lives in the native continuous P-2 space.
        // Extend F to redundant quadrature values by projecting its input;
        // otherwise off-space pressure modes are a spurious Newton nullspace.
        m_pressure->FwdTransLocalElmt(pressure,m_pressure->UpdateCoeffs());
        m_pressure->AvgAssemble(); m_pressure->GlobalToLocal();
        m_pressure->BwdTrans(m_pressure->GetCoeffs(),pressure);
        if (!m_hasOutflow)
        {
            const double mean = m_pressure->Integral(pressure)/m_pressure->Integral(ones);
            for (int q = 0; q < nq; ++q) pressure[q] -= mean;
        }
        Array<OneD,NekDouble> wall(wallMetric.size());
        for (int q = 0; q < wall.size(); ++q) wall[q] = values[(nf+1)*nq+q];
        m_wallStrip->SetIterationState(wall);
    };
    double change = 0, correction = 0, div = 0, weakDiv = 0;
    double globalPressureScale = 1;
    int iteration = 0;
    int newtonLinearIterations = 0,newtonBacktracks = 0;
    int nonpositiveLinearImages = 0;
    bool positiveLinearImage = true;
    auto iterationMap = [&](bool report) {
        RANSTrialGuard trial(m_nonlinearTrial,m_wallStrip.get(),!report);
        Closure(state, gradient, diffusion, gain, loss, time);
        UpdateLogWall(state);
        UpdateWall(state, predictor, gradient, pressure, time, alpha);
        UpdateOutflow(pressure, time);
        if (StartupWeight(time) > 0 && !m_startupLoad.size())
        {
            // The prepared parallel field has zero explicit advection. Cache
            // its complete native weak residual, including Robin/outlet loads.
            // Fading this fixed load is a continuation to the unforced model,
            // rather than overwriting a field during transport.
            m_startupLoad = Fields(nf);
            AuditStage(state, predictor, diffusion, gain, loss, pressure, time, alpha, -1);
        }
        if (report && m_auditStage && iteration == 0)
            AuditStage(state, predictor, diffusion, gain, loss, pressure, time, alpha, 0);
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
            if (m_startupLoad.size() && StartupWeight(time) > 0)
                Vmath::Svtvp(nc, -StartupWeight(time), m_startupLoad[i], 1,
                             rhs, 1, rhs, 1);
            if (m_terrainWall && i < dim)
            {
                Vmath::Vcopy(nc,rhs,1,velocityRhs[i],1);
                if (i == dim-1) m_linearIterations += m_matrixFree->SolveVelocity(
                    velocityRhs,diffusion,alpha,next,m_preconditionedNewton ? &state : nullptr);
            }
            else SolveField(i, rhs, diffusion[i], mass, next[i], alpha,
                            m_preconditionedNewton ? &state[i] : nullptr);
            if (report && (!m_terrainWall || i >= dim-1) && m_session->DefinesParameter("RANSDebug") && m_session->GetComm()->GetRank() == 0)
            {
                std::cout << "RANS iteration " << iteration << " field " << i
                          << " input " << *std::min_element(input[i].begin(), input[i].end())
                          << " gain " << *std::min_element(gain[i].begin(), gain[i].end())
                          << " result [" << *std::min_element(next[i].begin(), next[i].end())
                          << ',' << *std::max_element(next[i].begin(), next[i].end()) << "]" << std::endl;
            }
        }
        // A Newton residual x-F(x) is defined even when its linear Picard
        // image F(x) has a signed scalar component. No closure is evaluated
        // on that image. Physical trial inputs and accepted outputs retain
        // their strict positivity guards; Picard keeps its historical guard.
        positiveLinearImage = CheckState(next,!m_useNewton);
        if (!positiveLinearImage) ++nonpositiveLinearImages;
        Fields velocity(dim);
        for (int i = 0; i < dim; ++i) velocity[i] = next[i];
        correction = CorrectPressure(velocity, diffusion[0], pressure, alpha, div, weakDiv);
        if (report && m_session->DefinesParameter("RANSDebug") && m_session->GetComm()->GetRank() == 0)
            std::cout << "RANS pressure " << iteration << " correction " << correction
                      << " divergence " << div << " weak " << weakDiv
                      << " wall_constraint " << (m_wallStrip ? m_wallStrip->ConstraintResidual() : 0) << " range ["
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
            difference = Maximum(difference);
            change = std::max(change, difference / std::max(1.0, scale));
            if (report && m_session->DefinesParameter("RANSDebug") && m_session->GetComm()->GetRank() == 0)
                std::cout << "RANS state change " << iteration << " field " << i
                          << " absolute " << difference << std::endl;
            for (int q = 0; q < nq; ++q)
                state[i][q] += m_relaxation * (next[i][q] - state[i][q]);
        }
        double pressureScale = 1;
        for (double p : pressure) pressureScale = std::max(pressureScale, std::abs(p));
        globalPressureScale = Maximum(pressureScale);
    };
    auto converged = [&]() {
        return positiveLinearImage && change < m_tolerance && correction/globalPressureScale < m_pressureTolerance &&
            weakDiv < m_tolerance && (!m_wallStrip || m_wallStrip->ConstraintResidual() < m_tolerance);
    };
    int maximumIterations = m_maxIterations;
    if (m_useNewton) m_session->LoadParameter("RANSNewtonMaxIterations",maximumIterations,60);
    for (; iteration < maximumIterations; ++iteration)
    {
        if (m_useNewton && iteration == 1)
        {
            // Balance the primal/reaction off-diagonal blocks of x-F(x).
            // This scaling changes only the Krylov/line-search norm; all
            // physical state, continuity and wall criteria stay independent.
            const double scale = iterationScales[0]*m_wallStrip->IterationReactionScale()/std::sqrt(m_relaxation);
            RequireEverywhere(std::isfinite(scale) && scale > 0,"invalid Newton reaction scale");
            iterationScales[nf+1] = differenceScales[nf+1] = scale;
            for (int q = 0; q < wallMetric.size(); ++q) iterationMetric[(nf+1)*nq+q] = wallMetric[q]/(scale*scale);
        }
        const auto oldIteration = m_useNewton ? pack() : std::vector<double>{};
        iterationMap(true);
        if (converged()) break;
        if (m_useNewton && iteration == 0 && !positiveLinearImage)
        {
            // Keep the physical incoming fields for the first Newton guess.
            // The first sweep still initializes pressure and wall reaction.
            auto seed = pack();
            for (int i = 0; i < nf*nq; ++i) seed[i] = oldIteration[i];
            unpack(seed);
        }
        if (m_useNewton && iteration > 0)
        {
            // Newton solves x-F(x)=0, where F is one preconditioned Picard /
            // pressure / wall-reaction sweep. Native GMRES handles the full
            // nonsymmetric derivative coupling, including normal N samples.
            const auto mapped = pack();
            std::vector<double> residual(oldIteration.size());
            for (int i = 0; i < residual.size(); ++i) residual[i] = oldIteration[i]-mapped[i];
            const double norm2 = dotIteration(residual,residual);
            RequireEverywhere(std::isfinite(norm2),"nonfinite coupled Newton residual norm");
            if (!(norm2 > 0))
            {
                unpack(mapped); continue;
            }
            const int size = residual.size();
            Array<OneD,NekDouble> linearRhs(size),linearStep(size,0.0);
            for (int i = 0; i < size; ++i) linearRhs[i] = -std::sqrt(iterationMetric[i])*residual[i];
            unpack(oldIteration);
            Closure(state,gradient,diffusion,gain,loss,time);
            const auto baseBranches = m_cmuBranches;
            RANSLinearAction jacobian;
            jacobian.action = [&](LibUtilities::NekSysOperators::InArrayType &v,
                                  LibUtilities::NekSysOperators::OutArrayType &result,const bool &central) {
                double vn = 0;
                for (double value : v) vn += value*value;
                m_session->GetComm()->AllReduce(vn,LibUtilities::ReduceSum);
                if (vn == 0) { Vmath::Zero(size,result,1); return; }
                std::vector<double> direction(size),candidate(size);
                for (int i = 0; i < size; ++i) direction[i] = v[i]/std::sqrt(iterationMetric[i]);
                // Bound the *physical* perturbation, independent of tiny
                // GLL quadrature weights and mesh refinement. An L2-only
                // difference step can perturb a fine-grid vertex by O(1).
                double scaledMaximum = 0;
                for (int q = 0; q < nq; ++q)
                {
                    double velocity = 0;
                    for (int f = 0; f < dim; ++f) velocity = std::hypot(velocity,direction[f*nq+q]);
                    scaledMaximum = std::max(scaledMaximum,velocity/differenceScales[0]);
                    for (int f = dim; f <= nf; ++f)
                        scaledMaximum = std::max(scaledMaximum,std::abs(direction[f*nq+q])/differenceScales[f]);
                }
                for (int q = 0; q < wallMetric.size(); ++q)
                    scaledMaximum = std::max(scaledMaximum,std::abs(direction[(nf+1)*nq+q])/differenceScales[nf+1]);
                double step = m_newtonDifference/Maximum(scaledMaximum);
                for (int q = 0; q < nq; ++q)
                    for (int f : {m_k,m_epsilon})
                    {
                        const int i = f*nq+q;
                        const double minimum = f == m_k ? m_minimumK : m_minimumEpsilon;
                        if (direction[i] < 0 || (central && direction[i] != 0))
                            step = std::min(step,.05*(oldIteration[i]-minimum)/std::abs(direction[i]));
                    }
                m_session->GetComm()->AllReduce(step,LibUtilities::ReduceMin);
                bool valid = false;
                std::vector<double> plus;
                for (int attempt = 0; attempt < 8 && !valid; ++attempt)
                {
                    for (int i = 0; i < size; ++i) candidate[i] = oldIteration[i]+step*direction[i];
                    try
                    {
                        m_cmuBranches = baseBranches;
                        RANSFlagGuard branch(m_freezeCmuBranch,m_model == RANSModel::CLSHighRe);
                        unpack(candidate); iterationMap(false); plus = pack();
                        if (central)
                        {
                            for (int i = 0; i < size; ++i) candidate[i] = oldIteration[i]-step*direction[i];
                            unpack(candidate); iterationMap(false);
                        }
                        valid = true;
                    }
                    catch (const std::exception &) { step *= .1; }
                }
                RequireEverywhere(valid,"no admissible finite-difference Newton direction");
                const auto value = pack();
                for (int i = 0; i < size; ++i)
                    result[i] = std::sqrt(iterationMetric[i])*(central ?
                        direction[i]-(plus[i]-value[i])/(2*step) :
                        (candidate[i]-value[i]-residual[i])/step);
            };
            LibUtilities::NekSysOperators operators;
            operators.DefineNekSysLhsEval(&RANSLinearAction::Apply,&jacobian);
            LibUtilities::NekSysKey key;
            key.m_NekLinSysMaxIterations = m_newtonKrylovDimension;
            key.m_LinSysMaxStorage = m_newtonKrylovDimension;
            key.m_KrylovMaxHessMatBand = 0;
            key.m_NekLinSysTolerance = .01;
            // Freeze the active CLS max-invariant branch only for Jacobian
            // probes. Centered differences then approximate one linear
            // generalized Jacobian; physical/line-search residuals use max.
            key.m_GMRESCentralDifference = m_model == RANSModel::CLSHighRe;
            key.m_NekLinSysLeftPrecon = false;
            key.m_NekLinSysRightPrecon = false;
            auto linear = LibUtilities::GetNekLinSysIterFactory().CreateInstance(
                "GMRES",m_session,m_session->GetComm(),size,key);
            linear->SetSysOperators(operators); linear->SetFlagWarnings(false);
            linear->SetRhsMagnitude(norm2);
            const int linearIterations = linear->SolveSystem(size,linearRhs,linearStep,0);
            newtonLinearIterations += linearIterations;
            if (m_session->DefinesParameter("RANSDebug"))
            {
                Array<OneD,NekDouble> actual(size);
                jacobian.Apply(linearStep,actual,key.m_GMRESCentralDifference);
                double error = 0;
                for (int i = 0; i < size; ++i) error += (actual[i]-linearRhs[i])*(actual[i]-linearRhs[i]);
                m_session->GetComm()->AllReduce(error,LibUtilities::ReduceSum);
                if (m_session->GetComm()->GetRank() == 0)
                    std::cout << "RANS Newton " << iteration << " GMRES " << linearIterations
                              << " relative residual " << std::sqrt(error/norm2) << std::endl;
            }
            std::vector<double> direction(size),candidate(size);
            for (int i = 0; i < size; ++i) direction[i] = linearStep[i]/std::sqrt(iterationMetric[i]);
            double step = 1;
            for (int q = 0; q < nq; ++q)
                for (int f : {m_k,m_epsilon})
                {
                    const int i = f*nq+q;
                    const double minimum = f == m_k ? m_minimumK : m_minimumEpsilon;
                    if (direction[i] < 0) step = std::min(step,.95*(oldIteration[i]-minimum)/(-direction[i]));
                }
            for (double value : direction) if (!std::isfinite(value)) step = 0;
            m_session->GetComm()->AllReduce(step,LibUtilities::ReduceMin);
            bool accepted = false;
            for (int attempt = 0; attempt < 20 && step > 0 && !accepted; ++attempt)
            {
                for (int i = 0; i < size; ++i) candidate[i] = oldIteration[i]+step*direction[i];
                try
                {
                    unpack(candidate); iterationMap(false);
                    const auto value = pack();
                    std::vector<double> nextResidual(size);
                    for (int i = 0; i < size; ++i) nextResidual[i] = candidate[i]-value[i];
                    const double nextNorm2 = dotIteration(nextResidual,nextResidual);
                    accepted = std::isfinite(nextNorm2) &&
                        (nextNorm2 < (1-1e-4*step)*norm2 || converged());
                }
                catch (const std::exception &) { accepted = false; }
                if (!accepted) { step *= .5; ++newtonBacktracks; }
            }
            RequireEverywhere(accepted,"coupled Newton line search failed");
            if (converged()) break;
            unpack(candidate);
        }

    }
    RequireEverywhere(iteration < maximumIterations, "coupled implicit stage failed to converge");
    CheckState(next);
    // Return the actual Helmholtz solutions. They have the prescribed traces
    // and correspond to the last fixed-point residual, with no profile reset.
    for (int i = 0; i < nf; ++i)
        Vmath::Vcopy(nq, next[i], 1, output[i], 1);
    if (m_auditStage)
    {
        Closure(next, gradient, diffusion, gain, loss, time);
        UpdateWall(next, predictor, gradient, pressure, time, alpha);
        UpdateOutflow(pressure, time);
        AuditStage(next, predictor, diffusion, gain, loss, pressure, time, alpha, 1);
    }
    else if (m_terrainWall)
    {
        Closure(next,gradient,diffusion,gain,loss,time);
        UpdateWall(next,predictor,gradient,pressure,time,alpha);
    }
    if (m_wallStrip)
    {
        m_wallStrip->WriteDiagnostics(time, next, gain, loss, pressure);
        m_wallStrip->WritePointDiagnostics(time, 1);
    }
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
                      << maxNu << ',' << live << ",0," << m_linearIterations
                      << ',' << newtonLinearIterations << ',' << newtonBacktracks
                      << ',' << StartupWeight(time) << ',' << nonpositiveLinearImages << '\n' << std::flush;
}
} // namespace Nektar
