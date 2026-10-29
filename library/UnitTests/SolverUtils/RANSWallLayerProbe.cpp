///////////////////////////////////////////////////////////////////////////////
// File: RANSWallLayerProbe.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
//
// Milestone 1 research prototype. Native Nektar++ segment transforms, derivatives
// and weak products above H, with a separate control-volume strip below H.
// Milestone 2 evaluates the general CLS point kernel on simple-shear inputs.
// The wall reference quadrature stays independent. This is not an installed
// IncNavierStokesSolver equation system or a full nonlinear URANS timestep.
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Turbulence/RANSWallStrip.h>
#include <SolverUtils/Turbulence/RANSModels.h>
#include <LibUtilities/Communication/Comm.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LocalRegions/SegExp.h>
#include <SpatialDomains/PointGeom.h>
#include <SpatialDomains/SegGeom.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Nektar;
using namespace Nektar::SolverUtils;

namespace
{
const auto modelReference = RANSNeutralReference::CLS();
const double qStar = modelReference.Q();
const double kappa = modelReference.Kappa();
constexpr double sigmaE = 1.3;
using Vector = std::vector<double>;
using Matrix = std::vector<Vector>;

void Require(bool condition, const std::string &message)
{
    if (!condition)
    {
        throw std::runtime_error(message);
    }
}

double Maximum(const Vector &v)
{
    double result = 0.0;
    for (double value : v)
    {
        result = std::max(result, std::abs(value));
    }
    return result;
}

Vector Multiply(const Matrix &a, const Vector &x)
{
    Vector result(x.size(), 0.0);
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        result[i] = std::inner_product(a[i].begin(), a[i].end(), x.begin(), 0.0);
    }
    return result;
}

// Small dense prototype solve with partial pivoting; no production solver path
// is replaced. Matrices are assembled using native Nektar++ weak products.
Vector Solve(Matrix a, Vector b)
{
    const int n = b.size();
    for (int j = 0; j < n; ++j)
    {
        int pivot = j;
        for (int i = j + 1; i < n; ++i)
        {
            if (std::abs(a[i][j]) > std::abs(a[pivot][j]))
            {
                pivot = i;
            }
        }
        Require(std::abs(a[pivot][j]) > 1e-20, "Singular prototype matrix");
        std::swap(a[j], a[pivot]);
        std::swap(b[j], b[pivot]);
        for (int i = j + 1; i < n; ++i)
        {
            const double factor = a[i][j] / a[j][j];
            for (int k = j + 1; k < n; ++k)
            {
                a[i][k] -= factor * a[j][k];
            }
            b[i] -= factor * b[j];
        }
    }
    Vector x(n);
    for (int i = n - 1; i >= 0; --i)
    {
        double sum = b[i];
        for (int j = i + 1; j < n; ++j)
        {
            sum -= a[i][j] * x[j];
        }
        x[i] = sum / a[i][i];
    }
    return x;
}

double DualNorm(const Vector &residual, const Matrix &mass, int constrained = -1)
{
    Vector r;
    Matrix m;
    for (int i = 0; i < residual.size(); ++i)
    {
        if (i == constrained)
        {
            continue;
        }
        r.push_back(residual[i]);
        m.emplace_back();
        for (int j = 0; j < residual.size(); ++j)
        {
            if (j != constrained)
            {
                m.back().push_back(mass[i][j]);
            }
        }
    }
    auto inverse = Solve(m, r);
    const double squared = std::inner_product(r.begin(), r.end(), inverse.begin(), 0.0);
    Require(std::isfinite(squared) && squared >= -1e-25, "Invalid mass-dual norm");
    return std::sqrt(std::max(0.0, squared));
}

void Sum(const LibUtilities::CommSharedPtr &comm, Vector &values)
{
    Array<OneD, NekDouble> buffer(values.size());
    std::copy(values.begin(), values.end(), buffer.begin());
    comm->AllReduce(buffer, LibUtilities::ReduceSum);
    std::copy(buffer.begin(), buffer.end(), values.begin());
}

void Sum(const LibUtilities::CommSharedPtr &comm, Matrix &matrix)
{
    const int n = matrix.size();
    Vector values;
    values.reserve(n * n);
    for (const auto &row : matrix)
    {
        values.insert(values.end(), row.begin(), row.end());
    }
    Sum(comm, values);
    for (int i = 0; i < n; ++i)
    {
        std::copy_n(values.begin() + i * n, n, matrix[i].begin());
    }
}

struct ReferenceWall
{
    // Fixed inputs from the note; general wall kernels belong to milestone 2.
    double ustar = 0.4, nu = 1.5e-5, roughness = 0.1;
    RANSWallStripGeometry geometry{2.0, 1.0, 4.0};
    double K = qStar * ustar * ustar, t0 = ustar * ustar;
    double A = std::pow(ustar, 3) / kappa;
    double alpha = 0.2295, intercept, zv, slope, ze;
    double J0, J1, P0, epsilonMeanA, cutoff, matchingVelocity;

    ReferenceWall()
    {
        geometry.Validate();
        const double hStar = roughness * std::sqrt(K) / nu;
        const double x = hStar / 70.0;
        const double exponent = std::max(0.5 - 0.4 * std::pow(x, 0.7),
                                        1.0 - 0.79 * std::pow(x, -0.28));
        intercept = 10.7 * (1.0 - std::pow(x, exponent));
        zv = nu * intercept / std::sqrt(K);
        slope = alpha * std::sqrt(K);
        ze = 5.1 * nu / std::sqrt(K);
        J0 = Integrate([&](double z) { return 1.0 / D(z); });
        J1 = Integrate([&](double z) { return std::max(z - roughness, 0.0) / D(z); });
        P0 = t0 * t0 / geometry.height * Integrate([&](double z) {
            return NutA(z) / (D(z) * D(z));
        });
        epsilonMeanA = Integrate([&](double z) { return EpsilonA(z); }) / geometry.height;
        cutoff = geometry.height * std::exp(1.0 - geometry.height * P0 / A);
        Require(0.0 < cutoff && cutoff <= geometry.pointP,
                "Reference fails the note's cutoff admissibility condition");
        matchingVelocity = t0 * J0;
    }

    double NutA(double z) const { return slope * std::max(0.0, z - zv); }
    double D(double z) const { return nu + NutA(z); }
    double U(double z) const
    {
        Require(z > 0.0, "Outer logarithm evaluated at the wall");
        return matchingVelocity + ustar / kappa * std::log(z / geometry.height);
    }
    double Epsilon(double z) const
    {
        Require(z > 0.0, "Outer dissipation evaluated at the wall");
        return A / z;
    }
    double EpsilonA(double z) const
    {
        return z < ze ? 2.0 * nu * K / (ze * ze) : std::pow(K, 1.5) / (2.55 * z);
    }
    double EpsilonC(double z) const { return A / std::max(z, cutoff); }

    double Integrate(const std::function<double(double)> &function,
                     double upper = -1.0, Vector extra = {}) const
    {
        if (upper < 0.0)
        {
            upper = geometry.height;
        }
        Vector breaks{0.0, upper};
        extra.insert(extra.end(), {roughness, zv, ze});
        for (double value : extra)
        {
            if (value > 0.0 && value < upper)
            {
                breaks.push_back(value);
            }
        }
        for (double z = std::min(ze, nu / slope); z < upper; z *= 2.0)
        {
            breaks.push_back(z);
        }
        std::sort(breaks.begin(), breaks.end());
        breaks.erase(std::unique(breaks.begin(), breaks.end()), breaks.end());
        const LibUtilities::PointsKey key(40, LibUtilities::eGaussGaussLegendre);
        auto points = LibUtilities::PointsManager()[key];
        auto x = points->GetZ();
        auto w = points->GetW();
        double result = 0.0;
        for (int segment = 1; segment < breaks.size(); ++segment)
        {
            const double a = breaks[segment - 1], b = breaks[segment];
            for (int i = 0; i < x.size(); ++i)
            {
                result += 0.5 * (b - a) * w[i] *
                    function(0.5 * (a + b + (b - a) * x[i]));
            }
        }
        return result;
    }
};

void CheckStripRows(const ReferenceWall &wall)
{
    const auto &g = wall.geometry;
    RANSWallStripState state{wall.U(g.pointP), wall.K, wall.Epsilon(g.pointP)};
    RANSWallStripRates rates{0.0, 0.0};
    RANSWallStripTrace trace{wall.U(g.height), wall.U(g.pointN), wall.K, wall.K,
                            wall.Epsilon(g.height)};
    RANSWallStripExchange exchange{wall.t0, wall.t0, 0.0, 0.0, 0.0, 0.0,
                                  wall.P0, wall.P0, state.dissipation,
                                  trace.dissipationH};
    const auto residual = AssembleRANSWallStrip(g, state, rates, trace, exchange);
    Require(std::abs(residual.momentum) < 1e-13 && std::abs(residual.k) < 1e-13 &&
            std::abs(residual.dissipationP) < 1e-13 &&
            std::abs(residual.velocityMatch) < 1e-13 &&
            std::abs(residual.kMatch) < 1e-13 &&
            std::abs(residual.dissipationH) < 1e-13,
            "Reference wall-strip row or constraint is not balanced");
    Require(std::abs(g.EliminateVelocityP(trace.velocityH, trace.velocityN) -
                     state.velocity) < 1e-13 &&
            std::abs(g.EliminateKP(trace.kH, trace.kN) - state.k) < 1e-13,
            "Wall-state constraint elimination failed");

    // Non-equilibrium, independent row check. Summing strip and outer constant
    // test rows cancels the interface flux for arbitrary states, not just Q0.
    rates = {0.37, -0.11};
    exchange = {0.31, -0.07, 0.09, 0.0, -0.13, 0.04,
                0.6, 0.2, state.dissipation, trace.dissipationH};
    auto r = AssembleRANSWallStrip(g, state, rates, trace, exchange);
    Require(std::abs(r.momentum + exchange.momentumH -
            (g.height * (rates.velocity + exchange.meanMomentumForcing) +
             exchange.tractionWall)) < 1e-13, "Momentum interface counted twice");
    Require(std::abs(r.k + exchange.kH -
            (g.height * (rates.k + exchange.meanKAdvection) + exchange.kWall -
             g.height * (exchange.meanProduction - exchange.meanDissipation))) <
            1e-13, "TKE interface or source counted twice");
    state.dissipation *= 1.01;
    Require(std::abs(AssembleRANSWallStrip(g, state, rates, trace, exchange)
                         .dissipationP) > 1e-4,
            "Dissipation algebraic row does not detect a perturbation");
    Require(std::abs(RANSWallStripMassResidual(g, 0.3, 0.0, -0.15)) < 1e-13 &&
            std::abs(RANSWallStripMassResidual(g, 0.3, 0.0, 0.0)) > 0.1,
            "Integrated strip continuity failed");
    bool rejected = false;
    try
    {
        RANSWallStripGeometry{2.0, 2.0, 4.0}.Validate();
    }
    catch (const std::invalid_argument &)
    {
        rejected = true;
    }
    Require(rejected, "Degenerate strip geometry was accepted");
}

struct Record
{
    int degree, subdivisions, quadrature, dofs, sampleOwner;
    double momentum, tke, dissipation, normalMomentum, continuity;
    double uError, epsilonError, matching, stationaryError;
    double timestepDefect, timestepBalance, stationaryStepDefect;
    double noYap, zeroTopEpsilonFlux, wrongTraction, pressureGradient, crossflow;
};

Record Column(const ReferenceWall &wall, int degree, int subdivisions,
              int extraQuadrature, const LibUtilities::CommSharedPtr &comm)
{
    // Refinement leaves P, H and N fixed. In particular N=4 is a shared vertex,
    // owned by the element on its right, rather than a moving quadrature node.
    Vector boundaries;
    const Vector coarse{2.0, 4.0, 8.0, 16.0, 30.0};
    for (int e = 0; e + 1 < coarse.size(); ++e)
    {
        for (int part = 0; part < subdivisions; ++part)
        {
            boundaries.push_back(coarse[e] + (coarse[e + 1] - coarse[e]) *
                                 part / subdivisions);
        }
    }
    boundaries.push_back(coarse.back());
    const int elements = boundaries.size() - 1, modes = degree + 1;
    const int n = elements * degree + 1, nq = 2 * degree + extraQuadrature;
    Matrix mass(n, Vector(n, 0.0)), stiffness(n, Vector(n, 0.0));
    Vector ru(n, 0.0), rk(n, 0.0), re(n, 0.0), rn(n, 0.0), yapLoad(n, 0.0);
    Vector rc(n, 0.0), pressureLoad(n, 0.0), crossflowLoad(n, 0.0);
    Vector coefficients(n, 0.0), counts(n, 0.0), sample(n, 0.0);
    double errorU = 0.0, errorE = 0.0;
    int sampleOwner = -1;

    for (int e = 0; e < elements; ++e)
    {
        if (e % comm->GetSize() != comm->GetRank())
        {
            continue;
        }
        const double a = boundaries[e], b = boundaries[e + 1];
        SpatialDomains::PointGeom v0(1, 2 * e, a, 0.0, 0.0);
        SpatialDomains::PointGeom v1(1, 2 * e + 1, b, 0.0, 0.0);
        SpatialDomains::SegGeom geometry(e, 1, {&v0, &v1});
        LibUtilities::PointsKey points(nq, LibUtilities::eGaussLobattoLegendre);
        LibUtilities::BasisKey basis(LibUtilities::eModified_A, modes, points);
        LocalRegions::SegExp expansion(basis, &geometry);
        std::vector<int> map(modes);
        map[0] = e;
        map[1] = e + 1;
        for (int j = 2; j < modes; ++j)
        {
            map[j] = elements + 1 + e * (modes - 2) + j - 2;
        }
        Array<OneD, NekDouble> z(nq), u(nq), epsilon(nq), gradU(nq), gradE(nq);
        Array<OneD, NekDouble> cu(modes), ce(modes), work(nq), local(modes);
        expansion.GetCoords(z);
        for (int i = 0; i < nq; ++i)
        {
            u[i] = wall.U(z[i]);
            epsilon[i] = wall.Epsilon(z[i]);
        }
        expansion.FwdTransBndConstrained(u, cu);
        expansion.FwdTransBndConstrained(epsilon, ce);
        expansion.BwdTrans(cu, u);
        expansion.BwdTrans(ce, epsilon);
        expansion.PhysDeriv(0, u, gradU);
        expansion.PhysDeriv(0, epsilon, gradE);
        for (int j = 0; j < modes; ++j)
        {
            coefficients[map[j]] += cu[j];
            counts[map[j]] += 1.0;
        }
        if (a <= wall.geometry.pointN && wall.geometry.pointN < b)
        {
            sampleOwner = comm->GetRank();
            Array<OneD, NekDouble> xi(1);
            xi[0] = 2.0 * (wall.geometry.pointN - a) / (b - a) - 1.0;
            for (int j = 0; j < modes; ++j)
            {
                sample[map[j]] = expansion.PhysEvaluateBasis(xi, j);
            }
        }
        Array<OneD, NekDouble> fu(nq), fk(nq), fe(nq), se(nq), normal(nq), yy(nq);
        for (int i = 0; i < nq; ++i)
        {
            Require(epsilon[i] > 0.0, "Projected dissipation became nonpositive");
            RANSTensor gradient{};
            gradient[0][2]=gradU[i];
            const auto closure=EvaluateRANSClosure(
                RANSModel::CLSHighRe,wall.K,epsilon[i],gradient,z[i]);
            const double nut=closure.eddyViscosity, production=closure.production;
            const double Y=closure.yap;
            Require(std::isfinite(nut) && nut > 0.0 && std::isfinite(production) &&
                    std::isfinite(Y), "Nonfinite simple-shear constitutive state");
            fu[i] = nut * gradU[i];
            fk[i] = production - epsilon[i];
            fe[i] = nut / sigmaE * gradE[i];
            se[i] = closure.epsilonSource;
            // Nonlinear normal stress after absorbing 2k/3 into pressure.
            normal[i] = closure.quadraticStress[2][2]+closure.cubicStress[2][2];
            yy[i] = Y;
            errorU = std::max(errorU, std::abs(u[i] - wall.U(z[i])) / wall.ustar);
            errorE = std::max(errorE, std::abs(epsilon[i] / wall.Epsilon(z[i]) - 1.0));
        }
        auto assemble = [&](Vector &global, double scale) {
            for (int j = 0; j < modes; ++j)
            {
                global[map[j]] += scale * local[j];
            }
        };
        expansion.IProductWRTDerivBase(0, fu, local); assemble(ru, 1.0);
        expansion.IProductWRTBase(fk, local); assemble(rk, -1.0);
        expansion.IProductWRTDerivBase(0, fe, local); assemble(re, 1.0);
        expansion.IProductWRTBase(se, local); assemble(re, -1.0);
        expansion.IProductWRTDerivBase(0, normal, local); assemble(rn, -1.0);
        expansion.IProductWRTBase(yy, local); assemble(yapLoad, 1.0);

        // Actual weak continuity and normal-momentum controls. An affine W
        // departure has nonzero divergence; an affine reduced pressure creates
        // a normal force. Neither is hidden by the flat-reference zeros.
        Array<OneD, NekDouble> normalVelocity(nq, 0.0), divergence(nq);
        expansion.PhysDeriv(0, normalVelocity, divergence);
        expansion.IProductWRTBase(divergence, local); assemble(rc, 1.0);
        for (int i = 0; i < nq; ++i)
        {
            normalVelocity[i] = 0.001 * (z[i] - wall.geometry.height);
            work[i] = 0.01 * wall.t0 / wall.geometry.height *
                      (z[i] - wall.geometry.height);
        }
        expansion.PhysDeriv(0, normalVelocity, divergence);
        expansion.IProductWRTBase(divergence, local); assemble(crossflowLoad, 1.0);
        expansion.IProductWRTDerivBase(0, work, local); assemble(pressureLoad, -1.0);

        // Native mass and frozen reference diffusion matrices for the strip
        // elimination and backward-Euler consistency experiment.
        for (int j = 0; j < modes; ++j)
        {
            Array<OneD, NekDouble> unit(modes, 0.0), phi(nq), derivative(nq);
            unit[j] = 1.0;
            expansion.BwdTrans(unit, phi);
            expansion.IProductWRTBase(phi, local);
            for (int i = 0; i < modes; ++i)
            {
                mass[map[i]][map[j]] += local[i];
            }
            expansion.PhysDeriv(0, phi, derivative);
            for (int i = 0; i < nq; ++i)
            {
                work[i] = kappa * wall.ustar * z[i] * derivative[i];
            }
            expansion.IProductWRTDerivBase(0, work, local);
            for (int i = 0; i < modes; ++i)
            {
                stiffness[map[i]][map[j]] += local[i];
            }
        }
    }
    for (auto *v : {&ru, &rk, &re, &rn, &rc, &yapLoad, &pressureLoad,
                   &crossflowLoad, &coefficients, &counts, &sample})
    {
        Sum(comm, *v);
    }
    Sum(comm, mass);
    Sum(comm, stiffness);
    comm->AllReduce(sampleOwner, LibUtilities::ReduceMax);
    comm->AllReduce(errorU, LibUtilities::ReduceMax);
    comm->AllReduce(errorE, LibUtilities::ReduceMax);
    for (int i = 0; i < n; ++i)
    {
        coefficients[i] /= counts[i];
    }
    const int top = elements;
    ru[0] += wall.t0;
    // U at the top is Dirichlet. Its weak row is replaced, not counted as free.
    const double topFlux = -std::pow(wall.ustar, 4) / (sigmaE * boundaries.back());
    re[top] -= topFlux;
    const double qzz = -0.14666666666666666667 * wall.K;
    rn[0] -= qzz;
    rn[top] += qzz;
    const double rootLength = std::sqrt(boundaries.back() - boundaries.front());
    const double scaleU = wall.t0 / wall.geometry.height * rootLength;
    const double scaleK = std::pow(wall.ustar, 3) / wall.geometry.height * rootLength;
    const double scaleE = std::pow(wall.ustar, 4) /
                          (wall.geometry.height * wall.geometry.height) * rootLength;
    Record result{};
    result.degree = degree;
    result.subdivisions = subdivisions;
    result.quadrature = nq;
    result.dofs = n;
    result.sampleOwner = sampleOwner;
    result.momentum = DualNorm(ru, mass, top) / scaleU;
    result.tke = DualNorm(rk, mass) / scaleK;
    // The trace at H is the algebraic dissipation condition. No independent
    // epsilon Neumann condition or transport row is retained there.
    result.dissipation = DualNorm(re, mass, 0) / scaleE;
    result.normalMomentum = DualNorm(rn, mass) / scaleU;
    result.continuity = DualNorm(rc, mass) / (wall.ustar / wall.geometry.height * rootLength);
    result.uError = errorU;
    result.epsilonError = errorE;
    const double valueN = std::inner_product(sample.begin(), sample.end(),
                                            coefficients.begin(), 0.0);
    result.matching = std::abs(coefficients[0] - wall.geometry.VelocityAtH(
        wall.U(wall.geometry.pointP), valueN)) / wall.ustar;

    auto badE = re;
    for (int i = 0; i < n; ++i) { badE[i] += yapLoad[i]; }
    result.noYap = DualNorm(badE, mass, 0) / scaleE;
    badE = re;
    badE[top] += topFlux;
    result.zeroTopEpsilonFlux = DualNorm(badE, mass, 0) / scaleE;
    auto badU = ru;
    badU[0] += 0.01 * wall.t0;
    result.wrongTraction = DualNorm(badU, mass, top) / scaleU;
    pressureLoad[top] += 0.01 * wall.t0 / wall.geometry.height *
                         (boundaries.back() - wall.geometry.height);
    for (int i = 0; i < n; ++i) { pressureLoad[i] += rn[i]; }
    result.pressureGradient = DualNorm(pressureLoad, mass) / scaleU;
    result.crossflow = DualNorm(crossflowLoad, mass) /
                      (wall.ustar / wall.geometry.height * rootLength);

    // Exact elimination of UP and interface flux, preserving the strip's
    // inertia. c.U=(1-wN)UP. No ordinary CG cell occupies 0<z<H.
    const double weightN = wall.geometry.VelocityWeightN();
    const double weightP = 1.0 - weightN;
    Vector constraint(n, 0.0);
    constraint[0] = 1.0;
    for (int j = 0; j < n; ++j) { constraint[j] -= weightN * sample[j]; }
    Matrix coupledMass = mass, coupledStiffness = stiffness;
    for (int j = 0; j < n; ++j)
    {
        coupledMass[0][j] += wall.geometry.height / weightP * constraint[j];
    }
    coupledStiffness[0][0] += 1.0 / wall.J0;
    Matrix steadyMatrix = coupledStiffness;
    std::fill(steadyMatrix[top].begin(), steadyMatrix[top].end(), 0.0);
    steadyMatrix[top][top] = 1.0;
    Vector steadyRhs(n, 0.0);
    steadyRhs[top] = wall.U(boundaries.back());
    const Vector steady = Solve(steadyMatrix, steadyRhs);
    result.stationaryError = 0.0;
    for (int i = 0; i <= elements; ++i)
    {
        result.stationaryError = std::max(result.stationaryError,
            std::abs(steady[i] - wall.U(boundaries[i])) / wall.ustar);
    }
    const double dt = 0.125;
    Matrix stepMatrix = coupledStiffness;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j) { stepMatrix[i][j] += coupledMass[i][j] / dt; }
    }
    std::fill(stepMatrix[top].begin(), stepMatrix[top].end(), 0.0);
    stepMatrix[top][top] = 1.0;
    auto step = [&](const Vector &initial) {
        Vector rhs = Multiply(coupledMass, initial);
        for (auto &v : rhs) { v /= dt; }
        rhs[top] = wall.U(boundaries.back());
        return Solve(stepMatrix, rhs);
    };
    const Vector advanced = step(coefficients), steadyAdvanced = step(steady);
    Vector delta(n), steadyDelta(n), rate(n);
    for (int i = 0; i < n; ++i)
    {
        delta[i] = advanced[i] - coefficients[i];
        rate[i] = delta[i] / dt;
        steadyDelta[i] = steadyAdvanced[i] - steady[i];
    }
    result.timestepDefect = Maximum(delta) / wall.ustar;
    result.stationaryStepDefect = Maximum(steadyDelta) / wall.ustar;
    auto massRate = Multiply(mass, rate);
    const auto diffusion = Multiply(stiffness, advanced);
    const double topReaction = massRate[top] + diffusion[top];
    const double stripRate = std::inner_product(constraint.begin(), constraint.end(),
                                                rate.begin(), 0.0) / weightP;
    const double inventoryRate = std::accumulate(massRate.begin(),
                                                massRate.begin() + elements + 1, 0.0);
    result.timestepBalance = std::abs(inventoryRate + wall.geometry.height * stripRate +
                                   advanced[0] / wall.J0 - topReaction) / wall.t0;
    return result;
}

void WriteRecord(std::ostream &out, const Record &r)
{
    out << "{\"degree\":" << r.degree << ",\"subdivisions\":" << r.subdivisions
        << ",\"quadrature\":" << r.quadrature << ",\"dofs\":" << r.dofs
        << ",\"sample_owner\":" << r.sampleOwner
        << ",\"momentum\":" << r.momentum << ",\"k\":" << r.tke
        << ",\"epsilon\":" << r.dissipation
        << ",\"normal_momentum\":" << r.normalMomentum
        << ",\"continuity\":" << r.continuity
        << ",\"u_projection_error\":" << r.uError
        << ",\"epsilon_projection_error\":" << r.epsilonError
        << ",\"matching_constraint\":" << r.matching
        << ",\"linear_stationary_vertex_error\":" << r.stationaryError
        << ",\"linear_be_reference_defect\":" << r.timestepDefect
        << ",\"linear_be_conservation_error\":" << r.timestepBalance
        << ",\"linear_be_stationary_defect\":" << r.stationaryStepDefect
        << ",\"control_missing_yap\":" << r.noYap
        << ",\"control_zero_top_epsilon_flux\":" << r.zeroTopEpsilonFlux
        << ",\"control_wrong_wall_traction\":" << r.wrongTraction
        << ",\"control_pressure_gradient\":" << r.pressureGradient
        << ",\"control_crossflow_divergence\":" << r.crossflow << '}';
}
} // namespace

int main(int argc, char *argv[])
{
    LibUtilities::CommSharedPtr comm;
    try
    {
        bool mpi = false;
        std::string path = "rans-wall-layer.json";
        for (int i = 1; i < argc; ++i)
        {
            const std::string argument(argv[i]);
            if (argument == "--mpi") { mpi = true; }
            else if (argument == "--output" && i + 1 < argc) { path = argv[++i]; }
            else { throw std::invalid_argument("Use --mpi and/or --output PATH"); }
        }
        comm = LibUtilities::GetCommFactory().CreateInstance(
            mpi ? "ParallelMPI" : "Serial", argc, argv);
        Require(comm->GetSize() <= 16, "RANS project allocation exceeds 16 ranks");
        ReferenceWall wall;
        CheckStripRows(wall);
        const double mean = wall.Integrate([&](double z) {
            return wall.NutA(z) * std::pow(wall.t0 / wall.D(z), 2) - wall.EpsilonC(z);
        }, -1.0, {wall.cutoff}) / wall.geometry.height;
        const double moment = wall.Integrate([&](double z) {
            return (2.0 * z / wall.geometry.height - 1.0) *
                (wall.NutA(z) * std::pow(wall.t0 / wall.D(z), 2) - wall.EpsilonC(z));
        }, -1.0, {wall.cutoff}) / wall.geometry.height;
        Require(std::abs(mean) < 1e-11 && std::abs(moment) > 1e-3,
                "Wall mean and higher source moment experiment failed");
        const double auxiliaryP = wall.t0 * wall.Integrate([&](double z) {
            return 1.0 / wall.D(z);
        }, wall.geometry.pointP);
        Require(std::abs(auxiliaryP - wall.U(wall.geometry.pointP)) > 1e-3,
                "Auxiliary and resolved velocities were incorrectly identified");

        std::vector<Record> records;
        for (int degree : {2, 4, 6, 8, 10, 12, 14})
        {
            records.push_back(Column(wall, degree, 1, 8, comm));
        }
        for (int subdivisions : {2, 4, 8})
        {
            records.push_back(Column(wall, 4, subdivisions, 8, comm));
        }
        records.push_back(Column(wall, 10, 1, 20, comm));
        const auto &coarse = records.front();
        const auto &fine = records[6];
        Require(fine.momentum < 2e-6 && fine.tke < 2e-6 &&
                fine.dissipation < 2e-6 && fine.normalMomentum < 2e-6,
                "Fine analytical column residual failed");
        Require(fine.dissipation < coarse.dissipation * 1e-4 &&
                fine.momentum < coarse.momentum * 1e-4,
                "Analytical column does not converge with polynomial order");
        Require(fine.noYap > 100.0 * fine.dissipation &&
                fine.zeroTopEpsilonFlux > 100.0 * fine.dissipation &&
                fine.wrongTraction > 100.0 * fine.momentum &&
                fine.pressureGradient > 1e-3 && fine.crossflow > 1e-3 &&
                fine.continuity < 1e-13,
                "Residual controls failed to detect inconsistent equations or fluxes");
        for (const auto &r : records)
        {
            Require(r.matching < 1e-10 && r.timestepBalance < 1e-9 &&
                    r.stationaryStepDefect < 1e-9,
                    "Constraint, conservative timestep, or stationary-state test failed");
        }
        int failedOutput = 0;
        if (comm->GetRank() == 0)
        {
            std::ofstream output(path);
            output << std::setprecision(17)
                << "{\"schema\":1,\"prototype\":\"cv_strip_native_cg_column\",\"ranks\":"
                << comm->GetSize()
                << ",\"scope\":\"simple_shear_cls_outer_residual_and_frozen_linear_be\""
                << ",\"P\":" << wall.geometry.pointP << ",\"H\":" << wall.geometry.height
                << ",\"N\":" << wall.geometry.pointN << ",\"K\":" << wall.K
                << ",\"J0\":" << wall.J0 << ",\"P0\":" << wall.P0
                << ",\"original_epsilon_mean\":" << wall.epsilonMeanA
                << ",\"original_strip_k_residual\":" <<
                   wall.geometry.height * (wall.epsilonMeanA - wall.P0)
                << ",\"corrected_source_mean\":" << mean
                << ",\"corrected_source_linear_moment\":" << moment
                << ",\"epsilon_cutoff\":" << wall.cutoff
                << ",\"resolved_UP\":" << wall.U(wall.geometry.pointP)
                << ",\"auxiliary_UP\":" << auxiliaryP
                << ",\"strip_row_tests_passed\":true,\"records\":[\n";
            for (int i = 0; i < records.size(); ++i)
            {
                if (i) { output << ",\n"; }
                WriteRecord(output, records[i]);
            }
            output << "\n]}\n";
            failedOutput = !output.good();
            std::cout << std::setprecision(12)
                << "RANS wall-layer prototype passed: ranks=" << comm->GetSize()
                << " p14 residuals U=" << fine.momentum << " k=" << fine.tke
                << " eps=" << fine.dissipation << " normal=" << fine.normalMomentum
                << "; source mean=" << mean << " linear moment=" << moment << '\n';
        }
        comm->AllReduce(failedOutput, LibUtilities::ReduceMax);
        Require(!failedOutput, "Cannot write prototype report");
        comm->Finalise();
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "RANS wall-layer prototype: " << error.what() << '\n';
        if (comm) { comm->Finalise(); }
        return 1;
    }
}
