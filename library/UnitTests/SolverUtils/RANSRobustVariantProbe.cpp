///////////////////////////////////////////////////////////////////////////////
// File: RANSRobustVariantProbe.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <SolverUtils/Turbulence/RANSRobustVariant.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace Nektar::SolverUtils;

namespace
{
int checks = 0;

void Require(bool value, const std::string &message)
{
    ++checks;
    if (!value) { throw std::runtime_error(message); }
}

void Near(double a, double b, double tolerance, const std::string &message)
{
    Require(std::isfinite(a) && std::isfinite(b) &&
        std::abs(a-b) <= tolerance*std::max({1.,std::abs(a),std::abs(b)}),message);
}

std::string Quote(const std::string &value)
{
    std::string result = "\"";
    for (char c : value)
    {
        if (c == '"' || c == '\\') { result += '\\'; result += c; }
        else if (c == '\n') { result += "\\n"; }
        else if (c == '\r') { result += "\\r"; }
        else if (c == '\t') { result += "\\t"; }
        else { result += c; }
    }
    return result+'"';
}

void Vector(std::ostream &out, const std::vector<double> &values)
{
    out << '[';
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        if (i) { out << ','; }
        Require(std::isfinite(values[i]),"Nonfinite probe JSON value"); out << values[i];
    }
    out << ']';
}

void Tensor(std::ostream &out, const RANSTensor &value)
{
    out << '[';
    for (int i = 0; i < 3; ++i)
    {
        if (i) { out << ','; }
        Vector(out,{value[i][0],value[i][1],value[i][2]});
    }
    out << ']';
}

std::vector<double> Flatten(const RANSClosureValues &v)
{
    std::vector<double> result{v.cMu,v.eddyViscosity,v.strainInvariant,v.rotationInvariant,v.divergence,
        v.production,v.kSource,v.epsilonProduction,v.epsilonDestruction,v.yap,v.epsilonSource,
        v.momentumDiffusion,v.kDiffusion,v.epsilonDiffusion};
    for (const auto *tensor : {&v.reynoldsStress,&v.quadraticStress,&v.cubicStress})
        for (const auto &row : *tensor)for (double value : row) { result.push_back(value); }
    return result;
}

void Closure(std::ostream &out, const RANSClosureValues &v)
{
    out << "{\"stress\":"; Tensor(out,v.reynoldsStress);
    out << ",\"quadratic\":"; Tensor(out,v.quadraticStress);
    out << ",\"cubic\":"; Tensor(out,v.cubicStress);
    out << ",\"cMu\":" << v.cMu << ",\"nut\":" << v.eddyViscosity
        << ",\"strain\":" << v.strainInvariant << ",\"rotation\":" << v.rotationInvariant
        << ",\"divergence\":" << v.divergence << ",\"production\":" << v.production
        << ",\"k_source\":" << v.kSource << ",\"epsilon_production\":" << v.epsilonProduction
        << ",\"epsilon_destruction\":" << v.epsilonDestruction << ",\"yap\":" << v.yap
        << ",\"epsilon_source\":" << v.epsilonSource << ",\"momentum_diffusion\":" << v.momentumDiffusion
        << ",\"k_diffusion\":" << v.kDiffusion << ",\"epsilon_diffusion\":" << v.epsilonDiffusion << '}';
}

struct ClosureCase
{
    std::string name;
    RANSTensor gradient{}, direction{};
    double k = 1, epsilon = 1, distance = 1, nu = 0, dk = 0, de = 0;
    int branch = 0;
    RANSRobustVariant variant;
};

RANSClosureValues Evaluate(const ClosureCase &c, double step, bool original = false)
{
    auto G = c.gradient;
    for (int i = 0; i < 3; ++i)for (int j = 0; j < 3; ++j) { G[i][j] += step*c.direction[i][j]; }
    if (original)
    {
        return c.branch ? EvaluateRANSClosureLinearization(RANSModel::CLSHighRe,c.k+step*c.dk,
            c.epsilon+step*c.de,G,c.distance,{},c.nu,c.branch == 1) :
            EvaluateRANSClosure(RANSModel::CLSHighRe,c.k+step*c.dk,c.epsilon+step*c.de,
                G,c.distance,{},c.nu);
    }
    const RANSRobustVariant variant = c.variant;
    return c.branch ? EvaluateRANSRobustClosureLinearization(RANSModel::CLSHighRe,c.k+step*c.dk,
        c.epsilon+step*c.de,G,c.distance,{},c.nu,c.branch == 1,variant) :
        EvaluateRANSRobustClosure(RANSModel::CLSHighRe,c.k+step*c.dk,c.epsilon+step*c.de,
            G,c.distance,{},c.nu,variant);
}

void WriteClosure(std::ostream &out, const ClosureCase &c)
{
    const auto original = Evaluate(c,0,true), modified = Evaluate(c,0);
    long double production = 0;
    for (int i = 0; i < 3; ++i)for (int j = 0; j < 3; ++j)
        production -= static_cast<long double>(modified.reynoldsStress[i][j])*c.gradient[i][j];
    Near(modified.production,static_cast<double>(production),1e-13,"Production must use modified stress");
    Near(modified.kSource,modified.production-c.epsilon,1e-13,"Consistent modified k source");
    Near(modified.epsilonProduction,1.44*c.epsilon/c.k*modified.production,1e-13,"Consistent epsilon production");
    Near(modified.epsilonSource,modified.epsilonProduction-modified.epsilonDestruction+modified.yap,
        1e-13,"Consistent complete epsilon source");
    const double sink=(c.epsilon-std::min(0.,modified.production))/c.k;
    Near(std::max(0.,modified.production)-sink*c.k,modified.kSource,1e-13,"Signed Picard k source");
    const double epsilonRate=RANSParameters{}.cEpsilon2*c.epsilon/c.k;
    Near(modified.epsilonProduction+modified.yap-epsilonRate*c.epsilon,
        modified.epsilonSource,1e-13,"Picard epsilon destruction reaction");
    if (c.name == "expansion_negative_production")
    {
        Require(modified.production < 0 && sink > c.epsilon/c.k,"Negative production remains a sink");
    }
    if (c.variant.Original())
    {
        Require(Flatten(original) == Flatten(modified),"Original default must be exactly unchanged");
    }
    Require(original.cMu == modified.cMu && original.eddyViscosity == modified.eddyViscosity &&
        original.momentumDiffusion == modified.momentumDiffusion && original.kDiffusion == modified.kDiffusion &&
        original.epsilonDiffusion == modified.epsilonDiffusion && original.quadraticStress == modified.quadraticStress,
        "The robust wrapper must preserve the original coefficients and quadratic stress");
    if (c.name.find("reference_") == 0)
    {
        const auto before = Flatten(original), after = Flatten(modified);
        for (std::size_t i = 0; i < before.size(); ++i)
        { Near(before[i],after[i],1e-13,"Common beta and reference strain must preserve reference values"); }
    }
    out << "{\"name\":" << Quote(c.name) << ",\"branch\":"
        << Quote(c.branch == 0 ? "physical" : c.branch == 1 ? "strain" : "rotation")
        << ",\"beta\":" << c.variant.beta << ",\"gamma\":" << c.variant.gamma
        << ",\"reference_strain\":"; Tensor(out,c.variant.referenceStrain);
    out << ",\"k\":" << c.k << ",\"epsilon\":" << c.epsilon
        << ",\"distance\":" << c.distance << ",\"nu\":" << c.nu << ",\"gradient\":";
    Tensor(out,c.gradient); out << ",\"direction\":"; Tensor(out,c.direction);
    out << ",\"direction_k\":" << c.dk << ",\"direction_epsilon\":" << c.de
        << ",\"original\":"; Closure(out,original); out << ",\"modified\":"; Closure(out,modified);
    out << ",\"finite_differences\":[";
    bool first = true;
    for (double step : {1e-3,1e-4,1e-5,1e-6})
    {
        if (!first) { out << ','; } first = false;
        out << "{\"step\":" << step << ",\"plus\":"; Closure(out,Evaluate(c,step));
        out << ",\"minus\":"; Closure(out,Evaluate(c,-step)); out << '}';
    }
    out << "]}";
}

void Controls(std::ostream &out, const RANSCappedDestructionControls &c)
{
    out << "{\"theta\":" << c.theta << ",\"delta\":" << c.delta
        << ",\"smooth\":" << (c.smooth ? "true" : "false")
        << ",\"budget_relative_tolerance\":" << c.budgetRelativeTolerance << '}';
}

void Input(std::ostream &out, const RANSCappedDestructionInput &v)
{
    out << "{\"sampled_k\":" << v.sampledK << ",\"pre_e\":"; Vector(out,v.preE);
    out << ",\"pre_d\":"; Vector(out,v.preD); out << ",\"k\":"; Vector(out,v.k);
    out << ",\"metric_weights\":"; Vector(out,v.metricWeights); out << '}';
}

void Output(std::ostream &out, const RANSCappedDestructionValues &v)
{
    out << "{\"sampled_k\":" << v.sampledK << ",\"budget_e\":" << v.budgetE
        << ",\"budget_d\":" << v.budgetD << ",\"budget_mismatch\":" << v.budgetMismatch
        << ",\"b\":" << v.b << ",\"m\":" << v.m << ",\"lambda\":" << v.lambda
        << ",\"capped\":" << (v.capped ? "true" : "false")
        << ",\"dissipation_budget\":" << v.dissipationBudget
        << ",\"destruction_budget\":" << v.destructionBudget
        << ",\"redistribution_budget\":" << v.redistributionBudget;
    for (const auto &[name,values] : std::vector<std::pair<std::string,const std::vector<double>*>>{
        {"limiter",&v.limiter},{"e",&v.dissipation},{"d",&v.destruction},{"c",&v.redistribution},
        {"rate",&v.rate},{"rate_bound",&v.rateBound}})
    { out << ',' << Quote(name) << ':'; Vector(out,*values); }
    out << '}';
}

RANSCappedDestructionInput Shift(const RANSCappedDestructionInput &base,
    const RANSCappedDestructionDirection &direction, double step)
{
    auto result = base; result.sampledK += step*direction.sampledK;
    for (std::size_t i = 0; i < result.k.size(); ++i)
    {
        if (!direction.preE.empty()) { result.preE[i] += step*direction.preE[i]; }
        if (!direction.preD.empty()) { result.preD[i] += step*direction.preD[i]; }
        if (!direction.k.empty()) { result.k[i] += step*direction.k[i]; }
        if (!direction.metricWeights.empty()) { result.metricWeights[i] += step*direction.metricWeights[i]; }
    }
    return result;
}

struct SourceCase
{
    std::string name;
    RANSCappedDestructionInput input;
    RANSCappedDestructionControls controls;
    RANSCappedDestructionDirection direction;
    bool oneSided = false;
};

void WriteSource(std::ostream &out, const SourceCase &c)
{
    const auto value = EvaluateRANSCappedDestruction(c.input,c.controls);
    const auto derivative = RANSCappedDestructionDirectional(c.input,c.direction,c.controls);
    Near(value.destructionBudget,value.lambda*value.budgetD,1e-13,"Normalized destruction budget");
    Near(value.dissipationBudget,value.lambda*value.budgetE,1e-13,"Dependent dissipation budget");
    Near(value.redistributionBudget,value.lambda*value.budgetMismatch,1e-13,"Zero-budget redistribution to cutoff tolerance");
    for (std::size_t i = 0; i < value.rate.size(); ++i)
    {
        Require(value.rate[i] >= 0 && value.rate[i] <= value.rateBound[i]*(1+1e-13),"Bounded nonnegative destruction rate");
        const auto point = EvaluateRANSCappedDestructionPoint(c.input.preE[i],c.input.preD[i],
            c.input.k[i],c.input.sampledK,value,c.controls);
        Near(point.dissipation,value.dissipation[i],1e-14,"Point dependent profile");
        Near(point.destruction,value.destruction[i],1e-14,"Point destruction profile");
        Near(point.redistribution,value.redistribution[i],1e-14,"Point redistribution profile");
        Near(point.rate*c.input.k[i],point.destruction,1e-14,"Picard destruction reaction");
        if (c.input.k[i] == 0)
        {
            Require(value.destruction[i] == 0,"Zero-k destruction must vanish");
            const double limit=c.controls.smooth?0.:
                c.input.preD[i]/(value.m*c.controls.theta*c.input.sampledK);
            Near(point.rate,limit,1e-14,"Picard one-sided rate at zero k");
        }
    }
    if (std::all_of(value.limiter.begin(),value.limiter.end(),[](double w){return w == 1;}))
    {
        Require(value.b == 1 && value.m == 1 && value.lambda == 1,"Exact all-ray plateau normalization");
        Require(value.dissipation == c.input.preE && value.destruction == c.input.preD,"Exact all-ray plateau profiles");
    }
    out << "{\"name\":" << Quote(c.name) << ",\"controls\":"; Controls(out,c.controls);
    out << ",\"input\":"; Input(out,c.input); out << ",\"value\":"; Output(out,value);
    RANSCappedDestructionInput direction;
    direction.sampledK = c.direction.sampledK; direction.preE = c.direction.preE;
    direction.preD = c.direction.preD; direction.k = c.direction.k; direction.metricWeights = c.direction.metricWeights;
    out << ",\"direction\":"; Input(out,direction);
    out << ",\"directional\":"; Output(out,derivative);
    out << ",\"one_sided\":" << (c.oneSided ? "true" : "false") << ",\"finite_differences\":[";
    bool first = true;
    for (double step : {1e-3,1e-4,1e-5,1e-6})
    {
        if (!first) { out << ','; } first = false;
        out << "{\"step\":" << step << ",\"plus\":";
        Output(out,EvaluateRANSCappedDestruction(Shift(c.input,c.direction,step),c.controls));
        if (!c.oneSided)
        { out << ",\"minus\":"; Output(out,EvaluateRANSCappedDestruction(Shift(c.input,c.direction,-step),c.controls)); }
        out << '}';
    }
    out << "]}";
}

std::vector<SourceCase> Sources()
{
    RANSCappedDestructionInput base;
    base.sampledK = .5; base.preD = {.4,.9,1.3,.7,.2};
    base.preE = {1.1,.7,.9,.6,.5};
    const std::vector<double> r{.15,.5,1.,1.5,1.9}, dr{.2,.3,.5,.5,.5};
    double E = 0, D = 0;
    for (std::size_t i = 0; i < r.size(); ++i)
    {
        base.metricWeights.push_back(dr[i]*(1+.05*r[i])*(1+.02*r[i]));
        E += base.metricWeights[i]*base.preE[i]; D += base.metricWeights[i]*base.preD[i];
    }
    for (double &e : base.preE) { e *= D/E; }
    std::vector<SourceCase> result;
    for (bool smooth : {false,true})
    {
        const std::vector<std::pair<std::string,std::vector<double>>> profiles{
            {"plateau",{1.,.8,.3,.9,1.2}}, {"active_uncapped",{.02,.08,.3,.16,.4}},
            {"capped_depletion",{.001,.002,.003,.004,.005}},
            {"local_zero",{0.,.25,.3,.4,.5}}, {"whole_zero",{0.,0.,0.,0.,0.}}};
        for (const auto &[name,ratios] : profiles)
        {
            SourceCase c; c.name = name+(smooth ? "_smooth" : "_linear");
            c.input = base; c.controls.smooth = smooth;
            for (double ratio : ratios) { c.input.k.push_back(base.sampledK*ratio); }
            c.direction.sampledK = .03*base.sampledK;
            for (std::size_t i = 0; i < ratios.size(); ++i)
            {
                c.direction.preE.push_back(.04*base.preE[i]); c.direction.preD.push_back(.04*base.preD[i]);
                c.direction.k.push_back(.001*(i+1)*base.sampledK);
                c.direction.metricWeights.push_back(0);
            }
            // This nonuniform metric direction preserves the common budget
            // exactly together with the common pre-profile scale direction.
            c.direction.metricWeights[0] = .01;
            c.direction.metricWeights[1] = -.01*(base.preE[0]-base.preD[0])/(base.preE[1]-base.preD[1]);
            c.oneSided = name == "local_zero" || name == "whole_zero";
            result.push_back(c);
        }
        for (bool upper : {false,true})
        {
            SourceCase c; c.name = std::string("limiter_kink_")+(upper ? "up" : "down")+(smooth ? "_smooth" : "_linear");
            c.input.sampledK = 1; c.input.preE = {1,1}; c.input.preD = {1,1};
            c.input.metricWeights = {.5,.5}; c.controls.theta = .25; c.controls.delta = .5;
            c.controls.smooth = smooth; c.input.k = {.25,.5};
            c.direction.k = {upper ? .1 : -.1,0}; c.oneSided = true; result.push_back(c);
            c.name = std::string("cap_kink_")+(upper ? "up" : "down")+(smooth ? "_smooth" : "_linear");
            c.input.k = {.125,.125}; c.direction.k = {upper ? .1 : -.1,upper ? .1 : -.1};
            result.push_back(c);
        }
    }
    return result;
}

std::vector<std::pair<std::string,std::string>> Rejections()
{
    std::vector<std::pair<std::string,std::string>> result;
    auto reject = [&](const std::string &name, const std::function<void()> &function)
    {
        try { function(); }
        catch (const std::exception &error) { result.emplace_back(name,error.what()); ++checks; return; }
        throw std::runtime_error("Invalid robust-kernel input was accepted: "+name);
    };
    RANSCappedDestructionInput base{1,{1,1},{1,1},{1,1},{.5,.5}};
    for (const std::string name : {"negative_k","zero_K","negative_E","negative_D","zero_weight","bad_budget","bad_shape","nonfinite"})
    {
        auto bad = base;
        if (name == "negative_k") { bad.k[0] = -.1; }
        if (name == "zero_K") { bad.sampledK = 0; }
        if (name == "negative_E") { bad.preE[0] = -1; }
        if (name == "negative_D") { bad.preD[0] = -1; }
        if (name == "zero_weight") { bad.metricWeights[0] = 0; }
        if (name == "bad_budget") { bad.preE[0] = 2; }
        if (name == "bad_shape") { bad.k.pop_back(); }
        if (name == "nonfinite") { bad.k[0] = std::numeric_limits<double>::quiet_NaN(); }
        reject(name,[&](){ EvaluateRANSCappedDestruction(bad); });
    }
    for (const std::string name : {"theta_zero","theta_one","delta_zero","delta_one"})
    {
        RANSCappedDestructionControls bad;
        if (name == "theta_zero") { bad.theta = 0; }
        if (name == "theta_one") { bad.theta = 1; }
        if (name == "delta_zero") { bad.delta = 0; }
        if (name == "delta_one") { bad.delta = 1; }
        reject(name,[&](){ EvaluateRANSCappedDestruction(base,bad); });
    }
    for (const std::string name : {"negative_beta","negative_gamma","asymmetric_Sref","traceful_Sref"})
    {
        RANSRobustVariant bad; bad.gamma = .1;
        if (name == "negative_beta") { bad.beta = -.1; }
        if (name == "negative_gamma") { bad.gamma = -.1; }
        if (name == "asymmetric_Sref") { bad.referenceStrain[0][1] = 1; }
        if (name == "traceful_Sref") { bad.referenceStrain[0][0] = 1; }
        reject(name,[&](){ EvaluateRANSRobustClosure(RANSModel::CLSHighRe,1,1,{},1,{},0,bad); });
    }
    reject("inconsistent_point_denominator",[&](){ EvaluateRANSCappedDestructionPoint(1,1,1,1,.2,.3); });
    auto zero = base; zero.k[0] = 0; RANSCappedDestructionDirection outside; outside.k = {-1,0};
    reject("outward_zero_k_direction",[&](){ RANSCappedDestructionDirectional(zero,outside); });
    zero=base;zero.preD={0,2};outside={};outside.preD={-1,0};
    reject("outward_zero_D_direction",[&](){ RANSCappedDestructionDirectional(zero,outside); });
    return result;
}
} // namespace

int main(int argc, char **argv)
{
    std::string path;
    try
    {
        if (argc != 3 || std::string(argv[1]) != "--output")
        { throw std::invalid_argument("Use --output REPORT.json"); }
        path = argv[2];
        const auto reference = RANSNeutralReference::CLS();
        RANSTensor G{}, dx{}, dz{}, inplane{}, Sref{};
        G[0][2] = 1; dx[1][0] = 1; dz[1][2] = 1; Sref[0][2] = Sref[2][0] = 1;
        const RANSVector xi{.45,0,std::sqrt(1-.45*.45)}, a{-xi[2],0,.45};
        for (int i = 0; i < 3; ++i)for (int j = 0; j < 3; ++j) { inplane[i][j] = a[i]*xi[j]; }
        std::vector<ClosureCase> cases;
        for (int branch : {0,1,2})for (double beta : {1.,.4,.45})for (double gamma : {0.,.125})
        {
            for (const auto &[name,direction] : std::vector<std::pair<std::string,RANSTensor>>{
                {"reference_transverse_x",dx},{"reference_transverse_z",dz},{"reference_inplane",inplane}})
            {
                ClosureCase c; c.name = name; c.gradient = G; c.direction = direction;
                c.k = reference.Q(); c.branch = branch;
                c.variant.beta = beta; c.variant.gamma = gamma; c.variant.referenceStrain = Sref;
                cases.push_back(c);
            }
        }
        const RANSTensor general{{{.1,.4,.8},{-.2,-.3,.5},{.3,-.6,.2}}};
        ClosureCase expansion;
        expansion.name = "expansion_negative_production";
        expansion.gradient[0][0] = expansion.gradient[2][2] = .01;
        expansion.variant.beta = .4;
        cases.push_back(expansion);
        for (int branch : {0,1,2})for (double beta : {1.,.4,.45})
        {
            ClosureCase c; c.name = "general_3d"; c.gradient = general; c.direction = inplane;
            c.k = .7; c.epsilon = .6; c.dk = .031; c.de = -.027; c.distance = .3; c.nu = 1.5e-5;
            c.branch = branch; c.variant.beta = beta; c.variant.gamma = .125; c.variant.referenceStrain = Sref;
            cases.push_back(c);
            c.name = "off_reference_q5_shear"; c.gradient = G; c.direction = G;
            c.k = 5; c.epsilon = 1; c.dk = c.de = c.nu = c.variant.gamma = 0;
            cases.push_back(c);
        }
        std::ostringstream out; out << std::setprecision(17) << "{\"schema\":1,\"scope\":"
            << Quote("explicit variant kernels and derivatives; no global/coupled stability qualification")
            << ",\"reference\":{\"q\":" << reference.Q() << ",\"kappa\":" << reference.Kappa()
            << ",\"cMu\":" << reference.Cmu() << "},\"closure_cases\":[";
        bool first = true;
        for (const auto &c : cases)
        {
            if (!first) { out << ','; } first = false;
            WriteClosure(out,c);
        }
        out << "],\"source_cases\":["; first = true;
        for (const auto &c : Sources())
        { if (!first) { out << ','; } first = false; WriteSource(out,c); }
        out << "],\"rejections\":["; first = true;
        for (const auto &[name,message] : Rejections())
        {
            if (!first) { out << ','; } first = false;
            out << "{\"name\":" << Quote(name) << ",\"message\":" << Quote(message) << '}';
        }
        out << "],\"checks\":" << checks << ",\"success\":true}\n";
        std::ofstream file(path); if (!file) { throw std::runtime_error("Cannot open robust variant report"); }
        file << out.str(); file.close();
        if (!file) { throw std::runtime_error("Cannot write robust variant report"); }
        std::cout << "Robust kernel probe: " << checks << " checks\n"; return 0;
    }
    catch (const std::exception &error)
    {
        if (!path.empty())
        {
            std::ofstream out(path); out << "{\"schema\":1,\"success\":false,\"checks\":"
                << checks << ",\"failure\":" << Quote(error.what()) << "}\n";
        }
        std::cerr << error.what() << '\n'; return 1;
    }
}
