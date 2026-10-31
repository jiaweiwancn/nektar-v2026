///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <SolverUtils/Turbulence/RANSWallGeometry.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using namespace Nektar::SolverUtils;
namespace
{
int checks = 0, rejected = 0;
double error = 0;
void Near(double a, double b, double tolerance = 5e-11)
{
    ++checks;
    const double e = std::abs(a-b)/std::max({1.0,std::abs(a),std::abs(b)});
    error = std::max(error,e);
    if (!std::isfinite(e) || e > tolerance) throw std::runtime_error("Geometry equality failed");
}
template<class F> void Reject(F f)
{
    bool bad = false;
    try { f(); } catch (const std::exception &) { bad = true; }
    ++checks; ++rejected;
    if (!bad) throw std::runtime_error("Invalid geometry was accepted");
}
}
int main(int argc, char **argv)
{
    try
    {
        RANSWallGeometry flat, plane, hill;
        plane.terrain = RANSTerrain::Plane; plane.slope = .4; plane.origin = .2;
        hill.terrain = RANSTerrain::CosineHill;
        for (const auto &geometry : {flat,plane,hill})
            for (int i = 0; i <= 100; ++i)
            {
                const double x = -.3+.6*i/100;
                const auto f = geometry.At(x);
                Near(f.tangent[0]*f.normal[0]+f.tangent[2]*f.normal[2],0);
                Near(f.tangent[0]*f.tangent[0]+f.tangent[2]*f.tangent[2],1);
                for (double eta : {.002,.008,.016})
                {
                    const auto p = geometry.Point(x,eta,3);
                    const auto found = geometry.Closest(p[0],p[2]);
                    Near(found.distance,eta); Near(found.x,x); Near(p[1],3);
                }
                geometry.CheckStrip(f,.008,.016);
                Near(f.Volume(.008),.008-f.curvature*.008*.008/2);
            }
        Near(hill.At(0).curvature,-.04*std::acos(-1.0)*std::acos(-1.0)/(.02));
        Near(hill.At(.1).z,0); Near(hill.At(-.1).tangent[2],0);
        Reject([&] { hill.CheckStrip(hill.At(0),.02,.04); });
        Reject([&] { hill.CheckStrip(hill.At(-.099),.008,.1); });
        auto invalid = hill; invalid.hillHalfLength = 0;
        Reject([&] { invalid.Validate(); });
        Reject([&] { hill.Closest(std::nan(""),1); });
        AWFSpatialTransition ramp{true,-.4,-.2};
        Near(ramp.Weight(-1),1); Near(ramp.Weight(1),0); Near(ramp.Weight(-.3),.5);
        for (int i = 0; i <= 100; ++i)
        {
            const double r = i/100.0;
            Near(ramp.Weight(ramp.begin+r*(ramp.end-ramp.begin)),1-10*r*r*r+15*std::pow(r,4)-6*std::pow(r,5));
        }
        AWFSpatialTransition bad{true,1,0}; Reject([&] { bad.Validate(); });
        Near(AWFSpatialTransition{}.Weight(-5),1);
        // A coordinate ramp changes chi only. Both closures must still recover
        // exactly the original traction/production and all three epsilon data.
        for (const auto &model : {RANSNeutralReference::Standard(),RANSNeutralReference::CLS()})
        {
            AWFGeometry g{1.5e-5,.008,.002,.004};
            EquilibriumWallFunction wall(g,model.Q()*.3*.3,model);
            const AWFVector u{4,-.2}, forcing{.03,-.01};
            const auto original = wall.Original().Evaluate(u,forcing);
            const auto zero = wall.EvaluateFixedWeight(u,forcing,ramp.Weight(0),AWFFailurePolicy::RequirePreservation,false);
            Near(zero.values.traction[0],original.traction[0],0);
            Near(zero.values.meanProduction,original.meanProduction,0);
            Near(zero.values.meanDissipation,original.meanDissipation,0);
            Near(zero.values.dissipationP,original.dissipationP,0);
            Near(zero.values.dissipationH,original.dissipationH,0);
        }
        std::ofstream out(argc > 1 ? argv[1] : "rans-wall-geometry.json");
        out << std::setprecision(17) << "{\"checks\":" << checks << ",\"rejections\":" << rejected
            << ",\"maximum_geometry_error\":" << error << ",\"metric_cases\":[";
        bool first = true;
        for (double curvature : {-20.0,0.0,20.0})
            for (double roughness : {0.0,.002,.005})
                for (bool cls : {false,true})
                {
                    const auto model = cls ? RANSNeutralReference::CLS() : RANSNeutralReference::Standard();
                    const double k = model.Q()*.3*.3;
                    EquilibriumWallFunction wall({1.5e-5,.008,roughness,.004},k,model);
                    const AWFVector u{4,-.2}, forcing{.03,-.01};
                    for (double chi : {0.0,.37,1.0})
                    {
                        const auto &r = wall.Reference();
                        const auto v = IntegrateAWFMetric(wall.Original(),u,forcing,curvature,chi,r.numerator,r.cutoff);
                        if (!first) out << ',';
                        first = false;
                        out << "{\"curvature\":" << curvature << ",\"roughness\":" << roughness
                            << ",\"k\":" << k << ",\"chi\":" << chi << ",\"A\":" << r.numerator
                            << ",\"cutoff\":" << r.cutoff << ",\"P\":" << v.meanProduction
                            << ",\"D\":" << v.meanDissipation << ",\"shear_x\":" << v.integratedShear[0] << '}';
                    }
                }
        out << "]}\n";
        if (!out.good()) throw std::runtime_error("Cannot write geometry receipt");
        std::cout << "Geometry checks: " << checks << ", rejected inputs: " << rejected << '\n';
    }
    catch (const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}
