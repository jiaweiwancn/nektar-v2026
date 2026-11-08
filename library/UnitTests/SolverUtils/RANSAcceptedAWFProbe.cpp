///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <SolverUtils/Turbulence/EquilibriumAWFProfiles.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Nektar::SolverUtils;
namespace
{
int checks = 0;
void Check(bool ok, const std::string &message)
{
    ++checks;
    if (!ok) throw std::runtime_error(message);
}
std::string Quote(const std::string &s)
{
    std::string out = "\"";
    for (char c : s)
    {
        if (c == '"' || c == '\\') out += '\\';
        if (c == '\n') out += "\\n";
        else out += c;
    }
    return out+'"';
}
template<class F> void Reject(F f, const char *message)
{
    bool rejected = false;
    try { f(); } catch (const std::invalid_argument &) { rejected = true; }
    Check(rejected, message);
}
// P_m(2r/H-1), including its derivative with respect to physical r. These
// moments span every ordinary cubic scalar test mode, and two higher modes.
std::array<double,2> TestMode(int m, double r, double H)
{
    const double x = 2*r/H-1;
    double previous = 1, p = x, dprevious = 0, dp = 1;
    if (m == 0) return {1,0};
    for (int j = 2; j <= m; ++j)
    {
        const double next = ((2*j-1)*x*p-(j-1)*previous)/j;
        const double dn = ((2*j-1)*(p+x*dp)-(j-1)*dprevious)/j;
        previous = p; p = next; dprevious = dp; dp = dn;
    }
    return {p,2*dp/H};
}
void WriteVector(std::ostream &out, const AWFVector &a)
{ out << '[' << a[0] << ',' << a[1] << ']'; }

int Error(const std::string &path, const std::string &status, const std::exception &error, int code)
{
    std::ofstream out(path);
    out << "{\"status\":" << Quote(status) << ",\"reason\":" << Quote(error.what())
        << ",\"checks_completed\":" << checks << "}\n";
    std::cerr << error.what() << '\n';
    return out.good() ? code : 1;
}
}

int main(int argc, char **argv)
{
    std::string output = "accepted-awf-kernel.json";
    double H = 2, h = .1, nu = 1.5e-5, ustar = .4, suppliedK = -1;
    bool hasSuppliedK = false;
    AWFNormalMetric metric;
    AWFProfileOptions options;
    options.formulation = AWFProfileFormulation::CLSCalibratedBudgetAWF;
    try
    {
        for (int j = 1; j < argc; ++j)
        {
            const std::string name = argv[j];
            if (j+1 == argc) throw std::invalid_argument("Missing accepted-AWF probe argument");
            const std::string value = argv[++j];
            if (name == "--output") output = value;
            else if (name == "--height") H = std::stod(value);
            else if (name == "--roughness") h = std::stod(value);
            else if (name == "--nu") nu = std::stod(value);
            else if (name == "--ustar") ustar = std::stod(value);
            else if (name == "--sampled-k") { suppliedK = std::stod(value); hasSuppliedK = true; }
            else if (name == "--kappa1") metric.kappa1 = std::stod(value);
            else if (name == "--kappa2") metric.kappa2 = std::stod(value);
            else if (name == "--order") options.quadratureOrder = std::stoi(value);
            else throw std::invalid_argument("Unknown accepted-AWF probe argument: "+name);
        }
        const auto model = RANSNeutralReference::CLS();
        const auto wall = CalibratedAWFParameters(model);
        if (!(std::isfinite(ustar) && ustar > 0))
            throw std::invalid_argument("Accepted-AWF probe requires positive friction velocity");
        const double K = hasSuppliedK ? suppliedK : model.Q()*ustar*ustar;
        const AWFProfileGeometry geometry{nu,H,h,H,H};
        const EquilibriumAWFProfiles p(geometry,K,model,metric,options,wall);
        const auto &reference = p.Reference();
        const auto &original = p.Original();
        const auto &integrals = original.Integrals();
        const double a = wall.cMu*wall.cLength*std::sqrt(K);
        const double d = nu-a*original.ViscosityInterceptHeight();
        const bool flat = metric.kappa1 == 0 && metric.kappa2 == 0;
        const bool sufficientPositive = flat && original.ViscosityInterceptHeight() <= 0 && d >= 2*nu;
        const double energyLhs = reference.traction*reference.matchingVelocity;
        const double flatProduction = reference.traction*reference.traction*integrals.productionWeight;
        const double energyRhs = nu*reference.traction*reference.traction*integrals.j2+flatProduction;
        const double budgetScale = std::max(reference.productionBudget,reference.minimumBudget);
        Check(std::abs(wall.cMu*model.Q()*model.Q()-1) < 2e-14, "cMu calibration identity");
        Check(std::abs(a/(model.Kappa()*reference.frictionVelocity)-1) < 2e-14, "AWF viscosity-slope calibration identity");
        Check(std::abs(std::pow(K,1.5)/wall.cLength/reference.numerator-1) < 2e-14, "Dissipation-tail calibration identity");
        Check(std::abs(reference.budgetResidual) <= 2e-11*budgetScale, "Reference cutoff does not close its physical budget");
        Check(p.ReferenceK(H) == K && p.ReferenceKFlux(0) == 0, "Reference endpoint value/ground flux");
        Check(std::abs(p.ReferenceKFlux(H)) <= 2e-10*budgetScale, "Reference matching flux is nonzero");
        Check(std::abs(energyLhs-energyRhs) <= 2e-12*std::max(energyLhs,energyRhs), "Counterfactual flat mechanical-energy identity");
        Check(std::abs(p.Dissipation(H)*H/reference.numerator-1) < 2e-14, "Dependent epsilon endpoint");

        auto wrongGeometry = geometry; wrongGeometry.sampleHeight = H/2;
        Reject([&]{EquilibriumAWFProfiles q(wrongGeometry,K,model,metric,options,wall);}, "Nonendpoint K sampler accepted");
        wrongGeometry = geometry; wrongGeometry.tailStart = H/2;
        Reject([&]{EquilibriumAWFProfiles q(wrongGeometry,K,model,metric,options,wall);}, "Nonendpoint dissipation tail accepted");
        Reject([&]{EquilibriumAWFProfiles q(geometry,K,model,metric,options);}, "Omitted calibrated wall parameters accepted");
        Reject([&]{p.Reweighted(.5);}, "Accepted AWF acquired a chi reweighting");
        auto wrongOptions = options; wrongOptions.requestedWeight = .5;
        Reject([&]{EquilibriumAWFProfiles q(geometry,K,model,metric,wrongOptions,wall);}, "Accepted AWF acquired a chi switch");
        wrongOptions = options; wrongOptions.failurePolicy = AWFFailurePolicy::UseOriginal;
        Reject([&]{EquilibriumAWFProfiles q(geometry,K,model,metric,wrongOptions,wall);}, "Accepted AWF acquired a fallback");

        std::vector<double> samples{0,1e-6,1e-4,.001,.005,.01,.02,.05,.1,.25,.5,1,H};
        samples.erase(std::remove_if(samples.begin(),samples.end(),[&](double r){return r>H;}),samples.end());
        samples.push_back(reference.cutoff);
        std::sort(samples.begin(),samples.end());
        samples.erase(std::unique(samples.begin(),samples.end()),samples.end());
        const AWFVector matchingVelocity{reference.matchingVelocity+.15,-.11}, forcing{.023,-.017};
        const auto forced = original.Evaluate(matchingVelocity,forcing);
        const auto uground = original.Velocity(0,matchingVelocity,forcing);
        const auto uendpoint = original.Velocity(H,matchingVelocity,forcing);
        for (int j = 0; j < 2; ++j)
        {
            Check(uground[j] == 0, "Forced momentum no slip");
            Check(std::abs(uendpoint[j]-matchingVelocity[j]) <= 2e-13*std::max(1.,std::abs(matchingVelocity[j])), "Forced matching velocity");
            const double predicted = (matchingVelocity[j]-forcing[j]*integrals.j1)/integrals.j0;
            Check(std::abs(forced.traction[j]-predicted) <= 2e-13*std::max(1.,std::abs(predicted)), "Calibrated forced traction");
        }

        std::ofstream out(output); out << std::setprecision(17);
        Check(out.good(), "Cannot create accepted-AWF probe output");
        out << "{\"status\":\"ok\",\"formulation\":\"CLSCalibratedBudgetAWF-v1\","
            << "\"reference_scope\":\"stationary flat reference or frozen metric normal-column only; not terrain equilibrium\","
            << "\"order\":" << options.quadratureOrder << ",\"geometry\":{\"height\":" << H
            << ",\"roughness\":" << h << ",\"molecular_nu\":" << nu << ",\"sample_height\":" << H
            << ",\"tail_start\":" << H << ",\"kappa1\":" << metric.kappa1 << ",\"kappa2\":" << metric.kappa2
            << ",\"minimum_stretch\":" << metric.minimumStretch << "},\"model\":{\"q\":" << model.Q()
            << ",\"kappa\":" << model.Kappa() << ",\"sigma_k\":" << model.Parameters().sigmaK
            << "},\"wall\":{\"cMu\":" << wall.cMu << ",\"cLength\":" << wall.cLength
            << ",\"smoothIntercept\":" << wall.smoothIntercept << "},\"reference\":{\"K\":" << K
            << ",\"friction_velocity\":" << reference.frictionVelocity << ",\"traction\":" << reference.traction
            << ",\"matching_velocity\":" << reference.matchingVelocity << ",\"numerator\":" << reference.numerator
            << ",\"production_budget\":" << reference.productionBudget << ",\"minimum_budget\":" << reference.minimumBudget
            << ",\"lambda\":" << reference.productionBudget/reference.numerator << ",\"cutoff\":" << reference.cutoff
            << ",\"cutoff_iterations\":" << reference.cutoffIterations << ",\"budget_residual\":" << reference.budgetResidual
            << ",\"J0\":" << integrals.j0 << ",\"J1\":" << integrals.j1 << ",\"J2\":" << integrals.j2
            << ",\"production_weight\":" << integrals.productionWeight << ",\"a\":" << a
            << ",\"zv\":" << original.ViscosityInterceptHeight() << ",\"d\":" << d
            << ",\"aerodynamic_z0\":" << H*std::exp(-model.Kappa()*reference.matchingVelocity/reference.frictionVelocity)
            << ",\"ground_k\":" << p.ReferenceK(0) << ",\"rough_sufficient_positive\":" << (sufficientPositive?"true":"false")
            << ",\"flat_energy_lhs\":" << energyLhs << ",\"flat_energy_rhs\":" << energyRhs << "},\"samples\":[";
        double previousK = std::numeric_limits<double>::infinity();
        bool positive = true, nonincreasing = true;
        for (std::size_t j = 0; j < samples.size(); ++j)
        {
            if (j) out << ',';
            const double r = samples[j], kval = p.ReferenceK(r);
            positive = positive && kval > 0;
            nonincreasing = nonincreasing && kval <= previousK+2e-11*K; previousK = kval;
            const double p0 = p.ReferenceProduction(r); const auto source = p.AtProduction(r,p0);
            Check(source.redistribution == 0 && p.RedistributionFlux(r) == 0, "Prescribed redistribution is not exactly zero");
            Check(source.kSource == p0-source.dissipation || std::abs(source.kSource-(p0-source.dissipation)) <= 2e-15*std::max(p0,source.dissipation), "Prescribed sink missing from source");
            out << "{\"r\":" << r << ",\"k\":" << kval << ",\"k_gradient\":" << p.ReferenceKGradient(r)
                << ",\"Q\":" << p.ReferenceKFlux(r) << ",\"D\":" << p.ReferenceKDiffusivity(r)
                << ",\"nu_t\":" << original.EddyViscosity(r) << ",\"p0\":" << p0 << ",\"E\":" << source.dissipation
                << ",\"redistribution\":" << source.redistribution << ",\"source\":" << source.kSource
                << ",\"J\":" << metric.Jacobian(r) << ",\"U\":" << original.Velocity(r,{reference.matchingVelocity,0},{0,0})[0]
                << ",\"U_gradient\":" << original.VelocityGradient(r,{reference.matchingVelocity,0},{0,0})[0] << '}';
        }
        if (sufficientPositive) Check(positive && nonincreasing, "Sufficient rough-reference positivity/monotonicity failed");
        out << "],\"sampled_reference_positive\":" << (positive?"true":"false")
            << ",\"sampled_reference_nonincreasing\":" << (nonincreasing?"true":"false") << ",\"quadrature\":[";
        const auto &quadrature = p.Quadrature();
        for (std::size_t j = 0; j < quadrature.size(); ++j)
        { if (j) out << ','; const auto &q = quadrature[j]; out << "{\"r\":" << q.r << ",\"weight\":" << q.weight << ",\"J\":" << q.jacobian << '}'; }
        out << "],\"weak_basis\":\"shifted Legendre P_m(2r/H-1), m=0..5\",\"weak_modes\":[";
        for (int m = 0; m <= 5; ++m)
        {
            long double diffusion = 0, source = 0, forcedAccepted = 0, forcedEP = 0;
            for (const auto &q : quadrature)
            {
                const auto phi = TestMode(m,q.r,H);
                const double p0 = p.ReferenceProduction(q.r), ec = p.Dissipation(q.r);
                const double actual = p.Production(q.r,matchingVelocity,forcing);
                diffusion -= static_cast<long double>(q.weight)*p.ReferenceKFlux(q.r)*phi[1];
                source += static_cast<long double>(q.weight)*q.jacobian*(p0-ec)*phi[0];
                forcedAccepted += static_cast<long double>(q.weight)*q.jacobian*(actual-ec)*phi[0];
                forcedEP += static_cast<long double>(q.weight)*q.jacobian*(actual-p0)*phi[0];
            }
            Check(std::abs(diffusion+source) <= 5e-10*budgetScale, "Ordinary diffusion does not balance a reference weak mode");
            if (m) out << ',';
            out << "{\"degree\":" << m << ",\"diffusion\":" << static_cast<double>(diffusion)
                << ",\"source\":" << static_cast<double>(source) << ",\"residual\":" << static_cast<double>(diffusion+source)
                << ",\"forced_accepted_source\":" << static_cast<double>(forcedAccepted)
                << ",\"forced_calibrated_ep_source\":" << static_cast<double>(forcedEP) << '}';
        }
        out << "],\"cutoff_budget_samples\":[";
        for (int j = 0; j < 3; ++j)
        {
            const double c = reference.cutoff+(H-reference.cutoff)*j/2;
            if (j) out << ',';
            out << "{\"cutoff\":" << c << ",\"integral\":" << p.DissipationBudget(c) << '}';
        }
        out << "],\"forced\":{\"matching_velocity\":"; WriteVector(out,matchingVelocity);
        out << ",\"forcing\":"; WriteVector(out,forcing);
        out << ",\"traction\":"; WriteVector(out,forced.traction);
        out << ",\"mean_production\":" << forced.meanProduction << ",\"samples\":[";
        for (std::size_t j = 0; j < samples.size(); ++j)
        {
            if (j) out << ',';
            const double r = samples[j];
            out << "{\"r\":" << r << ",\"u\":"; WriteVector(out,original.Velocity(r,matchingVelocity,forcing));
            out << ",\"gradient\":"; WriteVector(out,original.VelocityGradient(r,matchingVelocity,forcing)); out << '}';
        }
        out << "]},\"checks_completed\":" << checks << "}\n";
        Check(out.good(), "Could not finish accepted-AWF probe output");
        std::cout << "Accepted AWF kernel: " << checks << " checks; " << output << '\n';
        return 0;
    }
    catch (const std::domain_error &error) { return Error(output,"inapplicable",error,2); }
    catch (const std::invalid_argument &error) { return Error(output,"inapplicable",error,2); }
    catch (const std::range_error &error) { return Error(output,"inapplicable",error,2); }
    catch (const std::exception &error) { return Error(output,"verification_failed",error,1); }
}
