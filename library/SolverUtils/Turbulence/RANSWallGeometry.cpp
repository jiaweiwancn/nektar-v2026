///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <SolverUtils/Turbulence/RANSWallGeometry.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace Nektar::SolverUtils
{
namespace
{
void Require(bool valid, const char *message)
{
    if (!valid) throw std::invalid_argument(message);
}
double CapIntegral(double H, double A, double cutoff, double curvature)
{
    Require(std::isfinite(A) && A > 0 && std::isfinite(cutoff) && cutoff > 0,
            "Invalid metric dissipation distribution");
    const double zero = cutoff < H ? A*(1+std::log(H/cutoff)) : A*H/cutoff;
    const double first = cutoff < H ? A*(H-cutoff/2) : A*H*H/(2*cutoff);
    return zero-curvature*first;
}
}

double RANSWallFrame::Jacobian(double eta) const { return 1-curvature*eta; }
double RANSWallFrame::Volume(double height) const
{
    return height*(1-curvature*height/2);
}

void RANSWallGeometry::Validate() const
{
    Require(terrain == RANSTerrain::Flat || terrain == RANSTerrain::Plane ||
            terrain == RANSTerrain::CosineHill, "Unknown RANS terrain geometry");
    Require(std::isfinite(origin) && std::isfinite(centre) && std::isfinite(slope) &&
            std::isfinite(hillHeight) && hillHeight >= 0 &&
            std::isfinite(hillHalfLength) && hillHalfLength > 0,
            "Invalid RANS terrain geometry parameters");
}

RANSWallFrame RANSWallGeometry::At(double x) const
{
    Validate(); Require(std::isfinite(x), "Nonfinite wall foot coordinate");
    double z = origin, first = 0, second = 0;
    if (terrain == RANSTerrain::Plane)
    {
        first = slope; z += slope*(x-centre);
    }
    else if (terrain == RANSTerrain::CosineHill && std::abs(x-centre) < hillHalfLength)
    {
        const double a = std::acos(-1.0)/hillHalfLength;
        const double angle = a*(x-centre);
        z += .5*hillHeight*(1+std::cos(angle));
        first = -.5*hillHeight*a*std::sin(angle);
        second = -.5*hillHeight*a*a*std::cos(angle);
    }
    // The paper's graph is C1 at its feet. Use the flat-side curvature at the
    // exact feet; one-sided interior geometry is verified at interior points.
    const double length = std::hypot(1.0, first);
    return {x,z,0,second/(length*length*length),
            {1/length,0,first/length},{-first/length,0,1/length}};
}

RANSWallFrame RANSWallGeometry::Closest(double x, double z) const
{
    Validate(); Require(std::isfinite(x) && std::isfinite(z), "Nonfinite terrain query");
    if (terrain != RANSTerrain::CosineHill || hillHeight == 0)
    {
        const double a = terrain == RANSTerrain::Plane ? slope : 0;
        auto f = At((x+a*(z-origin)+a*a*centre)/(1+a*a));
        f.distance = (x-f.x)*f.normal[0]+(z-f.z)*f.normal[2];
        return f;
    }
    const double left = centre-hillHalfLength, right = centre+hillHalfLength;
    std::vector<double> candidates{left,right,centre};
    if (x <= left || x >= right) candidates.push_back(x);
    auto residual = [&](double p) {
        const auto f = At(p);
        return (p-x)*f.tangent[0]+(f.z-z)*f.tangent[2];
    };
    // Search the whole compact hill, not only the vertical projection. There
    // can be multiple stationary distances above concave parts of a graph.
    constexpr int panels = 128;
    for (int i = 0; i < panels; ++i)
    {
        double a = left+(right-left)*i/panels, b = left+(right-left)*(i+1)/panels;
        double fa = residual(a), fb = residual(b);
        if (fa == 0) candidates.push_back(a);
        if (fb == 0) candidates.push_back(b);
        if ((fa < 0) == (fb < 0)) continue;
        for (int j = 0; j < 60; ++j)
        {
            const double m = .5*(a+b), fm = residual(m);
            if ((fa < 0) == (fm < 0)) { a = m; fa = fm; }
            else b = m;
        }
        candidates.push_back(.5*(a+b));
    }
    double best = std::numeric_limits<double>::max();
    RANSWallFrame result{};
    for (double p : candidates)
    {
        const auto f = At(p);
        const double square = (x-f.x)*(x-f.x)+(z-f.z)*(z-f.z);
        if (square < best) { result = f; best = square; }
    }
    const double signedNormal = (x-result.x)*result.normal[0]+(z-result.z)*result.normal[2];
    result.distance = std::copysign(std::sqrt(best), signedNormal);
    return result;
}

RANSVector RANSWallGeometry::Point(double x, double eta, double span) const
{
    Require(std::isfinite(eta) && std::isfinite(span), "Nonfinite normal coordinate");
    const auto f = At(x);
    return {f.x+eta*f.normal[0],span,f.z+eta*f.normal[2]};
}

void RANSWallGeometry::CheckStrip(const RANSWallFrame &f, double H, double N, double limit) const
{
    Require(std::isfinite(H) && H > 0 && std::isfinite(N) && N > H &&
            std::isfinite(limit) && limit > 0 && limit < 1,
            "Invalid normal strip dimensions/curvature limit");
    Require(std::abs(f.curvature)*H <= limit && f.Jacobian(N) > .1,
            "Terrain strip exceeds curvature limit or has a folded normal offset");
    for (double eta : {H,N})
    {
        const auto p = Point(f.x,eta);
        const auto c = Closest(p[0],p[2]);
        Require(std::abs(c.x-f.x) <= 1e-9*std::max(1.0,std::abs(f.x)) &&
                std::abs(c.distance-eta) <= 1e-9*std::max(1.0,eta),
                "Terrain strip normal ray does not have a unique matching foot");
    }
}

void AWFSpatialTransition::Validate() const
{
    Require(!enabled || (std::isfinite(begin) && std::isfinite(end) && end > begin),
            "AWF spatial transition requires finite begin < end");
}
double AWFSpatialTransition::Weight(double x) const
{
    Validate(); Require(std::isfinite(x), "Nonfinite AWF transition coordinate");
    if (!enabled || x <= begin) return 1;
    if (x >= end) return 0;
    const double r = (x-begin)/(end-begin), u = 1-r;
    return u*u*u*(1+3*r+6*r*r);
}

AWFMetricValues IntegrateAWFMetric(const AnalyticalWallFunction &wall,
    const AWFVector &u, const AWFVector &forcing, double curvature,
    double weight, double A, double cutoff)
{
    const auto &g = wall.Geometry(); const double H = g.height;
    Require(std::isfinite(curvature) && 1-curvature*H > 0 &&
            std::isfinite(weight) && weight >= 0 && weight <= 1,
            "Invalid AWF strip metric or correction weight");
    const auto base = wall.Evaluate(u,forcing);
    const double volume = H*(1-curvature*H/2);
    const double originalA = std::pow(wall.K(),1.5)/wall.Parameters().cLength;
    double dissipation = (1-weight)*CapIntegral(H,originalA,wall.DissipationCutoff(),curvature);
    if (weight > 0) dissipation += weight*CapIntegral(H,A,cutoff,curvature);
    double firstProduction = 0;
    if (curvature != 0)
    {
        constexpr double x[8] = {.09501250983763744,.2816035507792589,.4580167776572274,.6178762444026438,
            .755404408355003,.8656312023878318,.9445750230732326,.9894009349916499};
        constexpr double w[8] = {.1894506104550685,.1826034150449236,.1691565193950025,.1495959888165767,
            .1246289712555339,.0951585116824928,.0622535239386479,.027152459411754095};
        std::vector<double> breaks{0,H};
        for (double z : {g.roughness,wall.ViscosityInterceptHeight()})
            if (z > 0 && z < H) breaks.push_back(z);
        double scale = g.molecularNu/(wall.Parameters().cMu*wall.Parameters().cLength*std::sqrt(wall.K()));
        for (double z = scale; z < H; z *= 2) breaks.push_back(z);
        std::sort(breaks.begin(),breaks.end());
        auto production = [&](double z) {
            const auto du = wall.VelocityGradient(z,u,forcing);
            return z*wall.EddyViscosity(z)*(du[0]*du[0]+du[1]*du[1]);
        };
        for (size_t i = 1; i < breaks.size(); ++i)
        {
            const double mid = .5*(breaks[i]+breaks[i-1]), half = .5*(breaks[i]-breaks[i-1]);
            for (int j = 0; j < 8; ++j)
                firstProduction += half*w[j]*(production(mid-half*x[j])+production(mid+half*x[j]));
        }
    }
    const double length = std::max(0.0,H-g.roughness);
    return {(H*base.meanProduction-curvature*firstProduction)/volume,dissipation/volume,
        {H*base.traction[0]+.5*forcing[0]*length*length,
         H*base.traction[1]+.5*forcing[1]*length*length}};
}
}
