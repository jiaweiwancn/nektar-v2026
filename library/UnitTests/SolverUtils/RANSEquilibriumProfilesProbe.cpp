///////////////////////////////////////////////////////////////////////////////
// File: RANSEquilibriumProfilesProbe.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <SolverUtils/Turbulence/EquilibriumAWFProfiles.h>
#include <StdRegions/StdSegExp.h>
#include <StdRegions/StdQuadExp.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using namespace Nektar;
using namespace Nektar::SolverUtils;
namespace
{
int checks = 0, rejected = 0;
void Check(bool ok, const char *message)
{
    ++checks;
    if (!ok) throw std::runtime_error(message);
}
void Near(double a, double b, double scale = 1)
{
    Check(std::isfinite(a) && std::isfinite(b) &&
          std::abs(a-b) <= 2e-10*std::max(scale, std::abs(b)), "Native EP-AWF identity");
}
template<class F> void Reject(F function)
{
    bool bad = false;
    try { function(); } catch (const std::exception &) { bad = true; }
    Check(bad, "Invalid EP-AWF input accepted"); ++rejected;
}
void Vector(std::ostream &out, const std::vector<double> &v)
{
    out << '[';
    for (std::size_t i = 0; i < v.size(); ++i) { if (i) out << ','; out << v[i]; }
    out << ']';
}
void Point(std::ostream &out, const AWFProfileValues &v)
{
    Vector(out, {v.production, v.referenceProduction, v.dissipation, v.redistribution, v.kSource});
}
LibUtilities::BasisKey Basis(int degree, LibUtilities::BasisType type = LibUtilities::eModified_A)
{
    return {type, static_cast<std::size_t>(degree+1),
        {static_cast<std::size_t>(std::max(12, degree+3)), LibUtilities::eGaussLobattoLegendre}};
}
struct NativeBasis
{
    std::vector<std::vector<double>> phi, derivative, probes;
};
NativeBasis EvaluateBasis(StdRegions::StdSegExp &seg, const EquilibriumAWFProfiles &p)
{
    NativeBasis b;
    const double H = p.Geometry().height;
    for (int mode = 0; mode < seg.GetNcoeffs(); ++mode)
    {
        Array<OneD, NekDouble> coefficients(seg.GetNcoeffs(), 0.0), values(seg.GetTotPoints()),
            derivative(seg.GetTotPoints()), coord(1);
        coefficients[mode] = 1;
        seg.BwdTrans(coefficients, values);
        seg.PhysDeriv(0, values, derivative);
        b.phi.emplace_back(); b.derivative.emplace_back(); b.probes.emplace_back();
        for (const auto &q : p.Quadrature())
        {
            coord[0] = 2*q.r/H-1;
            b.phi.back().push_back(seg.PhysEvaluateBasis(coord, mode));
            b.derivative.back().push_back(2/H*seg.PhysEvaluate(coord, derivative));
        }
        for (double x : {-.81, .15, .77})
        {
            coord[0] = x;
            b.probes.back().push_back(seg.PhysEvaluateBasis(coord, mode));
            b.probes.back().push_back(2/H*seg.PhysEvaluate(coord, derivative));
        }
    }
    return b;
}
struct Case
{
    std::string name;
    AWFProfileGeometry geometry;
    AWFNormalMetric metric;
    double k;
    RANSNeutralReference model;
};
void WriteCase(std::ostream &out, const Case &c, int order, double chi)
{
    AWFProfileOptions options; options.quadratureOrder = order; options.requestedWeight = chi;
    EquilibriumAWFProfiles p(c.geometry, c.k, c.model, c.metric, options);
    const auto &r = p.Reference(); const auto &g = c.geometry;
    const AWFVector u{r.matchingVelocity, 0}, fa{-.02, 0}, fb{.03, .01};
    out << "{\"name\":\"" << c.name << "\",\"order\":" << order << ",\"chi\":" << chi
        << ",\"k\":" << c.k << ",\"q\":" << c.model.Q() << ",\"kappa\":" << c.model.Kappa()
        << ",\"geometry\":"; Vector(out, {g.molecularNu, g.height, g.roughness, g.sampleHeight, g.tailStart});
    out << ",\"metric\":"; Vector(out, {c.metric.kappa1, c.metric.kappa2});
    out << ",\"reference\":";
    Vector(out, {r.frictionVelocity, r.traction, r.matchingVelocity, r.numerator,
        r.productionBudget, r.minimumBudget, r.margin, r.reynolds, r.cutoff, r.weight,
        r.budgetResidual, static_cast<double>(r.cutoffIterations)});
    out << ",\"quadrature\":[";
    std::vector<double> source, meanOnly, correction, flux;
    for (std::size_t j = 0; j < p.Quadrature().size(); ++j)
    {
        const auto &q = p.Quadrature()[j];
        if (j) out << ',';
        Vector(out, {q.r, q.weight, q.jacobian});
        const auto v = p.Evaluate(q.r, u, {0, 0});
        source.push_back(v.kSource); meanOnly.push_back(v.production-v.dissipation);
        correction.push_back(v.redistribution); flux.push_back(p.RedistributionFlux(q.r));
        Check(q.weight > 0 && q.jacobian > 0, "Nonpositive split quadrature");
        Near(v.production-v.dissipation+v.redistribution, v.kSource,
             std::max(v.dissipation, std::abs(v.kSource)));
        if (chi == 1) Check(v.kSource == 0, "Reference source did not cancel before integration");
    }
    out << "],\"samples\":[";
    std::vector<double> heights{0, .01*g.height, .1*g.height, .5*g.height, g.height};
    if (r.cutoff > 0) for (double z : {.25*r.cutoff, r.cutoff, 1.1*r.cutoff})
        if (z <= g.height) heights.push_back(z);
    std::sort(heights.begin(), heights.end());
    for (std::size_t j = 0; j < heights.size(); ++j)
    {
        if (j) out << ',';
        const double z = heights[j];
        out << "{\"r\":" << z << ",\"zero\":"; Point(out, p.Evaluate(z, u, {0, 0}));
        out << ",\"forced\":"; Point(out, p.Evaluate(z, u, fb));
        out << ",\"supplied_negative_production\":"; Point(out, p.AtProduction(z, -.05));
        out << ",\"flux\":" << p.RedistributionFlux(z) << '}';
    }
    out << "],\"bases\":[";
    int count = 0;
    for (int degree : {2, 4, 6, 8})
    {
        StdRegions::StdSegExp seg(Basis(degree));
        const auto b = EvaluateBasis(seg, p);
        const auto loads = p.SourceLoads(b.phi, u, {0, 0});
        const auto red = p.RayLoads(b.phi, correction), old = p.RayLoads(b.phi, meanOnly);
        auto fl = p.RayLoads(b.derivative, flux);
        Array<OneD, NekDouble> coord(1);
        for (std::size_t i = 0; i < fl.size(); ++i)
        {
            coord[0] = 1;
            fl[i] = -fl[i]+seg.PhysEvaluateBasis(coord, i)*c.metric.Jacobian(g.height)*
                             p.RedistributionFlux(g.height);
            Near(red[i], fl[i], r.productionBudget);
            if (chi == 1) Check(loads[i] == 0, "Nonzero native reference source moment");
        }
        Array<OneD, NekDouble> nativeK(seg.GetTotPoints());
        const auto nodes = seg.GetBasis(0)->GetZ();
        for (int j = 0; j < nativeK.size(); ++j)
        {
            const double s = (nodes[j]+1)/2;
            nativeK[j] = c.k*(1+.2*s+.1*s*s);
        }
        const double sample = g.SampleK([&](double z) {
            coord[0] = 2*z/g.height-1;
            return seg.PhysEvaluate(coord, nativeK);
        });
        const double s = g.sampleHeight/g.height;
        Near(sample, c.k*(1+.2*s+.1*s*s), c.k);
        if (count++) out << ',';
        out << "{\"degree\":" << degree << ",\"sampled_polynomial_k\":" << sample
            << ",\"source\":"; Vector(out, loads);
        out << ",\"forced_a\":"; Vector(out, p.SourceLoads(b.phi, u, fa));
        out << ",\"forced_b\":"; Vector(out, p.SourceLoads(b.phi, u, fb));
        out << ",\"forced_scalar\":"; Vector(out, p.SourceLoads(b.phi, u, {.03, 0}));
        out << ",\"mean_only\":"; Vector(out, old);
        out << ",\"redistribution\":"; Vector(out, red);
        out << ",\"flux_action\":"; Vector(out, fl);
        out << ",\"basis_probes\":[";
        for (std::size_t j = 0; j < b.probes.size(); ++j)
        { if (j) out << ','; Vector(out, b.probes[j]); }
        out << "]}";
    }
    out << "],\"flux_endpoints\":"; Vector(out, {p.RedistributionFlux(0), p.RedistributionFlux(g.height)});
    // All active tensor-product basis rows at one physical wall ray. Tangential
    // integration and global interelement assembly are M3, not this ray action.
    if (c.name == "note" && order == 16 && chi == .37)
    {
        StdRegions::StdQuadExp quad(Basis(4), Basis(4));
        std::vector<std::vector<double>> phi(quad.GetNcoeffs());
        Array<OneD, NekDouble> x(2); x[0] = .27;
        for (int i = 0; i < quad.GetNcoeffs(); ++i)
            for (const auto &q : p.Quadrature())
            { x[1] = 2*q.r/g.height-1; phi[i].push_back(quad.PhysEvaluateBasis(x, i)); }
        out << ",\"quad_ray_loads\":"; Vector(out, p.SourceLoads(phi, u, fb));
    }
    if (c.name == "note" && chi == 1)
    {
        StdRegions::StdSegExp seg(Basis(6, LibUtilities::eOrtho_A));
        const auto b = EvaluateBasis(seg, p);
        auto loads = p.RayLoads(b.phi, meanOnly);
        out << ",\"orthonormal_mean_only\":"; Vector(out, loads);
        // Native Ortho_A is sqrt((2*j+1)/2)*L_j; the note tabulates L_j.
        for (std::size_t j = 0; j < loads.size(); ++j)
            loads[j] /= std::sqrt((2*j+1)/2.0);
        out << ",\"legendre_mean_only\":"; Vector(out, loads);
    }
    out << '}';
}
void Controls(std::ostream &out, double k)
{
    const auto model = RANSNeutralReference::CLS();
    const AWFProfileGeometry base{1.5e-5, 2, .1, 1, 1};
    for (auto g : {AWFProfileGeometry{0,2,.1,1,1}, {1.5e-5,2,2,1,1},
                   {1.5e-5,2,.1,0,1}, {1.5e-5,2,.1,3,1}, {1.5e-5,2,.1,1,3}})
        Reject([&] { EquilibriumAWFProfiles p(g, k, model); });
    Reject([&] { EquilibriumAWFProfiles p(base, 0, model); });
    Reject([&] { EquilibriumAWFProfiles p(base, k, model, {-.6, 0}); });
    AWFProfileOptions invalid; invalid.requestedWeight = 1.01;
    Reject([&] { EquilibriumAWFProfiles p(base, k, model, {}, invalid); });
    AWFProfileSwitch sw; sw.weights[6] = 0;
    Reject([&] { EquilibriumAWFProfileWeight(1, {}, 1, 1000, sw); });
    auto bad = base; bad.tailStart = 1e-12;
    Reject([&] { EquilibriumAWFProfiles p(bad, k, model); });
    AWFProfileOptions fallback; fallback.failurePolicy = AWFFailurePolicy::UseOriginal;
    EquilibriumAWFProfiles off(bad, k, model, {}, fallback);
    Check(off.Reference().weight == 0 && off.Reference().cutoffIterations == 0 &&
          off.Reference().reason == AWFFallbackReason::CutoffInadmissible, "Bad-budget fallback failed");
    fallback.failurePolicy = AWFFailurePolicy::RequirePreservation; fallback.requestedWeight = 0;
    EquilibriumAWFProfiles inactive(bad, k, model, {}, fallback);
    Check(inactive.Reference().cutoffIterations == 0, "Inactive inadmissible cutoff was solved");
    Near(inactive.Dissipation(.3), inactive.Original().OriginalDissipation(.3));
    auto low = base; low.height = 1e-4; low.roughness = 0;
    low.sampleHeight = low.tailStart = 5e-5;
    Reject([&] { EquilibriumAWFProfiles p(low, k, model); });
    fallback.requestedWeight = 1; fallback.failurePolicy = AWFFailurePolicy::UseOriginal;
    EquilibriumAWFProfiles lowFallback(low, k, model, {}, fallback);
    Check(lowFallback.Reference().reason == AWFFallbackReason::OutsideHighReRegime &&
          lowFallback.Reference().cutoffIterations == 0, "High-Re fallback failed");
    EquilibriumAWFProfiles p(base, k, model);
    Reject([&] { p.AtProduction(-1, 0); });
    Reject([&] { p.RayLoads({{1}}, {1}); });
    Reject([&] { base.SampleK([](double) { return -1; }); });
    out << "\"fallbacks\":[" << static_cast<int>(off.Reference().reason) << ','
        << off.Reference().cutoffIterations << ',' << inactive.Reference().cutoffIterations << ','
        << static_cast<int>(lowFallback.Reference().reason) << "],\"switch_samples\":[";
    const std::vector<std::array<double, 4>> inputs{{1,0,1,1000},{.4,0,1,1000},
        {1,.06,1,1000},{1,.2,1,1000},{1,0,.075,1000},{1,0,1,150},
        {0,0,-1,1},{1,0,-1,1000}};
    for (std::size_t j = 0; j < inputs.size(); ++j)
    {
        if (j) out << ',';
        auto x = inputs[j]; std::array<double,7> sensors{}; sensors[0] = x[1];
        out << '[' << x[0] << ',' << x[1] << ',' << x[2] << ',' << x[3] << ','
            << EquilibriumAWFProfileWeight(x[0], sensors, x[2], x[3]) << ']';
    }
    std::array<double,7> cancelled{}; cancelled[5] = 1; cancelled[6] = -1;
    Check(EquilibriumAWFProfileWeight(1, cancelled, 1, 1000) == 0,
          "Cancelling pressure and convection hid departures");
    out << "],\"cancelled_forcing_weight\":"
        << EquilibriumAWFProfileWeight(1, cancelled, 1, 1000);
    // A degree-five nonnegative polynomial with k'(0)=0 and a quartic zero
    // at r0 in the viscous sublayer, away from the physical K sample.
    auto g = base; g.roughness = 0;
    AWFProfileOptions zero; zero.requestedWeight = 0;
    EquilibriumAWFProfiles original(g, k, model, {}, zero);
    const double r0 = original.Original().ViscosityInterceptHeight()/3;
    StdRegions::StdSegExp seg(Basis(6));
    const auto nodes = seg.GetBasis(0)->GetZ();
    Array<OneD,NekDouble> values(seg.GetTotPoints()), first(seg.GetTotPoints()),
        second(seg.GetTotPoints()), coord(1);
    auto polynomial = [&](double r) { return k*std::pow(
        (r-r0)/(g.sampleHeight-r0), 4)*(1+4*r/r0)/(1+4*g.sampleHeight/r0); };
    for (int j = 0; j < values.size(); ++j) values[j] = polynomial(g.height*(nodes[j]+1)/2);
    seg.PhysDeriv(0, values, first); seg.PhysDeriv(0, first, second);
    const double sample = g.SampleK([&](double r) {
        coord[0] = 2*r/g.height-1; return seg.PhysEvaluate(coord, values);
    });
    coord[0] = 2*r0/g.height-1;
    const double value = seg.PhysEvaluate(coord, values);
    const double d1 = 2/g.height*seg.PhysEvaluate(coord, first);
    const double d2 = 4/(g.height*g.height)*seg.PhysEvaluate(coord, second);
    EquilibriumAWFProfiles sampledOriginal(g, sample, model, {}, zero);
    const auto v = sampledOriginal.Evaluate(r0, {sampledOriginal.Reference().matchingVelocity,0}, {0,0});
    coord[0] = -1;
    const double wallDerivative = 2/g.height*seg.PhysEvaluate(coord, first);
    const double rhs = v.kSource+g.molecularNu*d2, dt = 1e-6, next = value+dt*rhs;
    Near(sample,k,k); Near(value,0,k); Near(d1,0,k/g.height); Near(d2,0,k/(g.height*g.height));
    Near(wallDerivative,0,k/g.height);
    Check(original.Original().EddyViscosity(r0) == 0 && v.production == 0 &&
          rhs < 0 && next < 0, "Original sampled-K positivity failure not reproduced");
    out << ",\"positivity_counterexample\":{\"r0\":" << r0 << ",\"rK\":" << g.sampleHeight
        << ",\"sampled_k\":" << sample << ",\"polynomial_value\":" << value
        << ",\"polynomial_first\":" << d1 << ",\"polynomial_second\":" << d2
        << ",\"analytic_value\":0,\"analytic_first\":0,\"analytic_second\":0,\"wall_derivative\":0"
        << ",\"native_wall_derivative\":" << wallDerivative
        << ",\"production\":" << v.production << ",\"dissipation\":" << v.dissipation
        << ",\"source\":" << v.kSource << ",\"rhs\":" << rhs << ",\"dt\":" << dt
        << ",\"local_euler_value\":" << next << ",\"full_timestep\":false}";
}
}
int main(int argc, char **argv)
{
    try
    {
        if (argc != 3 || std::string(argv[1]) != "--output")
            throw std::invalid_argument("Use --output PATH");
        const auto cls = RANSNeutralReference::CLS();
        const double k = cls.Q()*.16;
        const AWFProfileGeometry g{1.5e-5,2,.1,1,1};
        std::vector<Case> cases{{"note",g,{},k,cls}};
        auto smooth = g; smooth.roughness = 0;
        cases.push_back({"smooth",smooth,{},k,cls});
        auto partial = g; partial.roughness = 10*g.molecularNu/std::sqrt(k);
        cases.push_back({"partial",partial,{},k,cls});
        cases.push_back({"curved_note",g,{.05,0},k,cls});
        cases.push_back({"metric_two_curvatures",g,{.05,.08},k,cls});
        cases.push_back({"metric_mixed",g,{.07,-.04},k,cls});
        cases.push_back({"small_K",g,{},.25*k,cls});
        cases.push_back({"large_K",g,{},4*k,cls});
        auto tail = g; tail.tailStart = tail.height;
        cases.push_back({"tail_at_H",tail,{},k,cls});
        auto sample = g; sample.sampleHeight = sample.height;
        cases.push_back({"sample_at_H",sample,{},k,cls});
        std::ofstream out(argv[2]); out << std::setprecision(17) << "{\"schema\":1,\"cases\":[";
        int count = 0;
        for (const auto &c : cases)
            for (int order : {16,32}) for (double chi : {0.,.37,1.})
            { if (count++) out << ','; WriteCase(out,c,order,chi); }
        const auto standard = RANSNeutralReference::Standard();
        for (int order : {16,32})
        {
            out << ',';
            WriteCase(out,{"standard_control",g,{},standard.Q()*.16,standard},order,1);
            ++count;
        }
        out << "],"; Controls(out,k);
        out << ",\"summary\":{\"checks\":" << checks << ",\"rejected_inputs\":" << rejected
            << ",\"cases\":" << count << "}}\n";
        out.close();
        Check(static_cast<bool>(out), "Cannot finish profile report");
        std::cout << "M2 native profiles: " << count << " cases, " << checks
                  << " checks; original inner-k positivity gate fails as specified.\n";
    }
    catch (const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
    return 0;
}
