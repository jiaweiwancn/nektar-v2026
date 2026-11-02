///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <SolverUtils/Turbulence/EquilibriumAWFTransition.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Nektar::SolverUtils
{
namespace
{
void Finite(double x)
{
    if (!std::isfinite(x)) throw std::invalid_argument("Nonfinite EP-AWF transition input");
}
void Positive(double x)
{
    Finite(x);
    if (!(x > 0)) throw std::invalid_argument("EP-AWF transition requires positive input");
}
double Cast(long double x)
{
    const double result = static_cast<double>(x);
    if (!std::isfinite(result)) throw std::range_error("Unrepresentable EP-AWF transition result");
    return result;
}
struct Jet
{
    long double v, d, dd;
    Jet(long double value = 0, long double first = 0, long double second = 0)
        : v(value), d(first), dd(second) {}
};
Jet Make(double v, double d = 0, double dd = 0)
{
    Finite(v); Finite(d); Finite(dd);
    return {v, d, dd};
}
Jet operator+(const Jet &a, const Jet &b) { return {a.v+b.v, a.d+b.d, a.dd+b.dd}; }
Jet operator-(const Jet &a, const Jet &b) { return {a.v-b.v, a.d-b.d, a.dd-b.dd}; }
Jet operator*(const Jet &a, const Jet &b)
{
    return {a.v*b.v, a.d*b.v+a.v*b.d, a.dd*b.v+2*a.d*b.d+a.v*b.dd};
}
Jet operator/(const Jet &a, const Jet &b)
{
    if (b.v == 0) throw std::domain_error("Zero EP-AWF transition denominator");
    const Jet reciprocal{1/b.v, -b.d/(b.v*b.v),
        (2*b.d*b.d-b.v*b.dd)/(b.v*b.v*b.v)};
    return a*reciprocal;
}
Jet Sqrt(const Jet &a)
{
    if (!(a.v > 0)) throw std::domain_error("Nonpositive EP-AWF transition scale");
    const long double value = std::sqrt(a.v);
    return {value, a.d/(2*value), a.dd/(2*value)-a.d*a.d/(4*value*value*value)};
}
Jet Gate(const Jet &x, const Jet &lo, const Jet &hi)
{
    if (!(lo.v < hi.v)) throw std::invalid_argument("Unordered EP-AWF transition gate");
    if (x.v <= lo.v) return {1};
    if (x.v >= hi.v) return {0};
    const Jet s = (x-lo)/(hi-lo), u = Jet(1)-s;
    return u*u*u*(Jet(1)+Jet(3)*s+Jet(6)*s*s);
}
using JVector = std::array<Jet, 3>;
JVector Vector(const RANSVector &v, const RANSVector &d)
{
    JVector result;
    for (int a = 0; a < 3; ++a) result[a] = Make(v[a], d[a]);
    return result;
}
JVector Subtract(const JVector &a, const JVector &b)
{
    JVector result;
    for (int j = 0; j < 3; ++j) result[j] = a[j]-b[j];
    return result;
}
Jet Dot(const JVector &a, const JVector &b)
{
    Jet result;
    for (int j = 0; j < 3; ++j) result = result+a[j]*b[j];
    return result;
}
JVector Tangent(const JVector &value, const JVector &normal)
{
    JVector result;
    const Jet projection = Dot(value, normal);
    for (int j = 0; j < 3; ++j) result[j] = value[j]-normal[j]*projection;
    return result;
}
void Store(const Jet &j, double &v, double &d, double &dd)
{
    v = Cast(j.v); d = Cast(j.d); dd = Cast(j.dd);
}
void Store(const JVector &j, RANSVector &v, RANSVector &d, RANSVector &dd)
{
    for (int a = 0; a < 3; ++a) Store(j[a], v[a], d[a], dd[a]);
}
Jet SensorSquare(double x, double dx) { const Jet a = Make(x, dx); return a*a; }
Jet SensorSquare(const RANSVector &x, const RANSVector &dx)
{
    const JVector a = Vector(x, dx);
    return Dot(a, a);
}
void CheckWall(const RANSWallNormalMap &wall, const RANSVector &chart)
{
    for (double q : chart) Finite(q);
    if (wall.kind != RANSWallNormalMap::Kind::Flat &&
        wall.kind != RANSWallNormalMap::Kind::Plane &&
        wall.kind != RANSWallNormalMap::Kind::CosineHill)
        throw std::invalid_argument("EP-AWF transition requires a declared flat upstream support");
    // At validates the map parameters and dimension without relying on rho.
    wall.At({chart[0], chart[1], 0});
}
std::pair<double, double> Slope(const RANSWallNormalMap &wall, double x)
{
    if (wall.kind == RANSWallNormalMap::Kind::Flat ||
        wall.kind == RANSWallNormalMap::Kind::Plane) return {wall.slopeX, 0};
    if (wall.kind != RANSWallNormalMap::Kind::CosineHill || std::abs(x) >= wall.lengthX)
        return {0, 0};
    const double wave = std::acos(-1.0)/wall.lengthX;
    return {-.5*wall.amplitude*wave*std::sin(wave*x),
            -.5*wall.amplitude*wave*wave*std::cos(wave*x)};
}
std::vector<std::pair<long double, long double>> Gauss(int order)
{
    std::vector<std::pair<long double, long double>> rule;
    for (int j = 0; j < order/2; ++j)
    {
        long double z = std::cos(std::acos(-1.0L)*(j+.75L)/(order+.5L)), derivative = 0;
        for (int it = 0; it < 64; ++it)
        {
            long double p0 = 1, p1 = z;
            for (int n = 2; n <= order; ++n)
            {
                const long double p = ((2*n-1)*z*p1-(n-1)*p0)/n;
                p0 = p1; p1 = p;
            }
            derivative = order*(z*p1-p0)/(z*z-1);
            const long double step = p1/derivative;
            z -= step;
            if (std::abs(step) <= 8*std::numeric_limits<long double>::epsilon()) break;
            if (it == 63) throw std::runtime_error("EP-AWF arc quadrature did not converge");
        }
        const long double weight = 2/((1-z*z)*derivative*derivative);
        for (int side : {-1, 1}) rule.emplace_back(side*z, weight);
    }
    return rule;
}
long double ArcPanel(const RANSWallNormalMap &wall, double a, double b, int depth = 0)
{
    static const auto rule16 = Gauss(16), rule32 = Gauss(32);
    const auto integrate = [&](const auto &rule)
    {
        long double sum = 0;
        for (const auto &q : rule)
        {
            const double x = Cast((a+b)/2+static_cast<long double>(b-a)*q.first/2);
            sum += q.second*std::hypot(1.0, Slope(wall, x).first);
        }
        return sum*(b-a)/2;
    };
    const long double coarse = integrate(rule16), fine = integrate(rule32);
    if (std::abs(fine-coarse) <= 2e-14L*std::max(1.0L, std::abs(fine))) return fine;
    if (depth == 20) throw std::runtime_error("EP-AWF arc quadrature did not resolve the wall curve");
    const double midpoint = a+(b-a)/2;
    return ArcPanel(wall, a, midpoint, depth+1)+ArcPanel(wall, midpoint, b, depth+1);
}
double ArcIntegral(const RANSWallNormalMap &wall, double upper)
{
    if (wall.kind != RANSWallNormalMap::Kind::CosineHill)
    {
        const double slope = Slope(wall, upper).first;
        return Cast(static_cast<long double>(upper)*std::hypot(1.0, slope));
    }
    // The only nonconstant segment is [-L,L]; evenness lets each signed
    // integral use [0,min(abs(q1),L)] plus its exactly flat remainder.
    const double sign = upper < 0 ? -1 : 1;
    const double end = std::abs(upper), inside = std::min(end, wall.lengthX);
    if (inside == 0) return 0;
    return Cast(sign*(ArcPanel(wall, 0, inside)+end-inside));
}
} // namespace

void AWFTransitionControls::Validate() const
{
    state.Validate(); Finite(start); Finite(end);
    if (!(start < end) || !std::isfinite(end-start) ||
        !std::isfinite(state.off*state.off) || !(state.on*state.on > 0) ||
        !(state.on*state.on < state.off*state.off))
        throw std::invalid_argument("Invalid EP-AWF physical transition interval or state thresholds");
}

AWFTransitionResult EvaluateEquilibriumAWFTransition(const AWFTransitionInput &i,
    const AWFTransitionControls &c, const AWFTransitionInput &d,
    const AWFTransitionControlDirection &dc)
{
    c.Validate();
    if (i.reynolds < 0) throw std::invalid_argument("Negative EP-AWF transition Reynolds number");
    const auto &s = i.sensors, &ds = d.sensors;
    const std::array<Jet, 7> squares{
        SensorSquare(s.velocity, ds.velocity), SensorSquare(s.velocityDifference, ds.velocityDifference),
        SensorSquare(s.normal, ds.normal), SensorSquare(s.k, ds.k), SensorSquare(s.epsilon, ds.epsilon),
        SensorSquare(s.pressure, ds.pressure), SensorSquare(s.convection, ds.convection)};
    Jet departure;
    for (int j = 0; j < 7; ++j)
        departure = departure+Make(c.state.weights[j], dc.weights[j])*squares[j];
    const Jet on = Make(c.state.on, dc.on), off = Make(c.state.off, dc.off);
    const Jet spatial = Gate(Make(i.coordinate, d.coordinate), Make(c.start, dc.start), Make(c.end, dc.end));
    const Jet state = Gate(departure, on*on, off*off);
    const Jet margin = Jet(1)-Gate(Make(i.budgetMargin, d.budgetMargin),
        Make(c.state.marginOff, dc.marginOff), Make(c.state.marginOn, dc.marginOn));
    const Jet reynolds = Jet(1)-Gate(Make(i.reynolds, d.reynolds),
        Make(c.state.reynoldsOff, dc.reynoldsOff), Make(c.state.reynoldsOn, dc.reynoldsOn));
    const Jet admissibility = margin*reynolds, chi = spatial*admissibility*state;
    AWFTransitionResult result;
    Store(chi, result.value.chi, result.derivative.chi, result.secondDerivative.chi);
    Store(spatial, result.value.spatial, result.derivative.spatial, result.secondDerivative.spatial);
    Store(state, result.value.state, result.derivative.state, result.secondDerivative.state);
    Store(admissibility, result.value.admissibility, result.derivative.admissibility, result.secondDerivative.admissibility);
    Store(margin, result.value.marginGate, result.derivative.marginGate, result.secondDerivative.marginGate);
    Store(reynolds, result.value.reynoldsGate, result.derivative.reynoldsGate, result.secondDerivative.reynoldsGate);
    Store(departure, result.value.d2, result.derivative.d2, result.secondDerivative.d2);
    return result;
}

AWFTransitionSensorResult EvaluateEquilibriumAWFTransitionSensors(
    const AWFTransitionPhysicalInput &i, const AWFTransitionPhysicalInput &d)
{
    Positive(i.sampledK); Positive(i.height); Positive(i.molecularNu); Positive(i.Q);
    Positive(i.kA); Positive(i.epsilonA); Positive(i.epsilonAReference);
    const JVector normal = Vector(i.normal, d.normal);
    if (std::abs(Dot(normal, normal).v-1) > 2e-12L)
        throw std::invalid_argument("EP-AWF transition requires a unit geometric normal");
    const Jet K = Make(i.sampledK, d.sampledK), H = Make(i.height, d.height);
    const Jet nu = Make(i.molecularNu, d.molecularNu), Q = Make(i.Q, d.Q);
    const Jet viscous = nu/H, scale2 = K/Q+viscous*viscous, scale = Sqrt(scale2);
    const JVector matching = Tangent(Subtract(Vector(i.matchingVelocity, d.matchingVelocity),
        Vector(i.matchingReference, d.matchingReference)), normal);
    const JVector deltaA = Subtract(Vector(i.velocityA, d.velocityA), Vector(i.velocityAReference, d.velocityAReference));
    const JVector deltaB = Subtract(Vector(i.velocityB, d.velocityB), Vector(i.velocityBReference, d.velocityBReference));
    const JVector difference = Tangent(Subtract(deltaB, deltaA), normal);
    const JVector pressure = Subtract(Vector(i.pressure, d.pressure), Vector(i.pressureReference, d.pressureReference));
    const JVector convection = Subtract(Vector(i.convection, d.convection), Vector(i.convectionReference, d.convectionReference));
    AWFTransitionSensorResult result;
    JVector u, du, p, a;
    for (int j = 0; j < 3; ++j)
    {
        u[j] = matching[j]/scale; du[j] = difference[j]/scale;
        p[j] = H*pressure[j]/scale2; a[j] = H*convection[j]/scale2;
    }
    Store(u, result.value.velocity, result.derivative.velocity, result.secondDerivative.velocity);
    Store(du, result.value.velocityDifference, result.derivative.velocityDifference, result.secondDerivative.velocityDifference);
    Store(p, result.value.pressure, result.derivative.pressure, result.secondDerivative.pressure);
    Store(a, result.value.convection, result.derivative.convection, result.secondDerivative.convection);
    Store(Dot(normal, deltaA)/scale, result.value.normal, result.derivative.normal, result.secondDerivative.normal);
    Store((Make(i.kA, d.kA)-K)/K, result.value.k, result.derivative.k, result.secondDerivative.k);
    const Jet epsilonReference = Make(i.epsilonAReference, d.epsilonAReference);
    Store((Make(i.epsilonA, d.epsilonA)-epsilonReference)/epsilonReference,
        result.value.epsilon, result.derivative.epsilon, result.secondDerivative.epsilon);
    return result;
}

AWFTransitionCoordinateResult EquilibriumAWFWallCoordinateDifferential(
    const RANSWallNormalMap &wall, const RANSVector &chart, const RANSVector &direction)
{
    CheckWall(wall, chart);
    for (double q : direction) Finite(q);
    const auto slope = Slope(wall, chart[0]);
    AWFTransitionCoordinateResult result;
    result.coordinate = ArcIntegral(wall, chart[0]);
    result.chartGradient[0] = std::hypot(1.0, slope.first);
    result.derivative = Cast(static_cast<long double>(result.chartGradient[0])*direction[0]);
    result.secondDerivative = Cast(static_cast<long double>(slope.first)*slope.second/
        result.chartGradient[0]*direction[0]*direction[0]);
    return result;
}

double EquilibriumAWFWallCoordinate(const RANSWallNormalMap &wall, const RANSVector &chart)
{
    return EquilibriumAWFWallCoordinateDifferential(wall, chart).coordinate;
}

void ValidateEquilibriumAWFTransitionSupport(const RANSWallNormalMap &wall,
    const AWFTransitionControls &controls)
{
    controls.Validate(); CheckWall(wall, {0, 0, 0});
    if (wall.kind == RANSWallNormalMap::Kind::CosineHill &&
        !(controls.end < EquilibriumAWFWallCoordinate(wall, {-wall.lengthX, 0, 0})))
        throw std::invalid_argument("EP-AWF spatial transition must finish before the incoming hill foot");
}

AWFTransitionCoordinateResult EquilibriumAWFAffineWallCoordinate(
    const RANSWallNormalMap &wall, const RANSVector &chart, const RANSVector &origin,
    const RANSVector &covector, const RANSVector &direction)
{
    CheckWall(wall, chart);
    if (wall.kind == RANSWallNormalMap::Kind::CosineHill)
        throw std::invalid_argument("Affine EP-AWF coordinate requires a globally flat/plane support");
    const auto ground = wall.At({chart[0], chart[1], 0});
    long double norm2 = 0, normal = 0, coordinate = 0;
    AWFTransitionCoordinateResult result;
    for (int a = 0; a < 3; ++a)
    {
        Finite(origin[a]); Finite(covector[a]); Finite(direction[a]);
        norm2 += static_cast<long double>(covector[a])*covector[a];
        normal += static_cast<long double>(covector[a])*ground.normal[a];
        coordinate += static_cast<long double>(covector[a])*(ground.point[a]-origin[a]);
        for (int q = 0; q < 2; ++q)
            result.chartGradient[q] += covector[a]*ground.jacobian[a][q];
    }
    if (!(norm2 > 0) || std::abs(normal) > 2e-12L*std::sqrt(norm2))
        throw std::invalid_argument("EP-AWF affine coordinate requires a nonzero tangent covector");
    result.coordinate = Cast(coordinate);
    result.derivative = Cast(static_cast<long double>(result.chartGradient[0])*direction[0]+
        static_cast<long double>(result.chartGradient[1])*direction[1]);
    return result;
}

AWFTransitionProfileResult BlendEquilibriumAWFTransition(const AWFTransitionResult &t,
    const AWFTransitionProfileInput &i, const AWFTransitionProfileInput &d,
    const AWFTransitionProfileInput &dd)
{
    const Jet chi = Make(t.value.chi, t.derivative.chi, t.secondDerivative.chi);
    if (!(chi.v >= 0 && chi.v <= 1)) throw std::invalid_argument("Invalid EP-AWF blend weight");
    const Jet original = Make(i.originalDissipation, d.originalDissipation, dd.originalDissipation);
    const Jet production = Make(i.production, d.production, dd.production);
    Jet dissipation = original, redistribution, source = production-original;
    if (chi.v != 0 || chi.d != 0 || chi.dd != 0)
    {
        const Jet corrected = Make(i.correctedDissipation, d.correctedDissipation, dd.correctedDissipation);
        const Jet reference = Make(i.referenceProduction, d.referenceProduction, dd.referenceProduction);
        dissipation = original+chi*(corrected-original);
        redistribution = chi*(corrected-reference);
        source = source+chi*(original-reference);
    }
    AWFTransitionProfileResult result;
    Store(dissipation, result.value.dissipation, result.derivative.dissipation, result.secondDerivative.dissipation);
    Store(redistribution, result.value.redistribution, result.derivative.redistribution, result.secondDerivative.redistribution);
    Store(source, result.value.kSource, result.derivative.kSource, result.secondDerivative.kSource);
    return result;
}
} // namespace Nektar::SolverUtils
