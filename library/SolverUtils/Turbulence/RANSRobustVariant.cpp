///////////////////////////////////////////////////////////////////////////////
// File: RANSRobustVariant.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <SolverUtils/Turbulence/RANSRobustVariant.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace Nektar::SolverUtils
{
namespace
{
using Real = long double;

void Require(bool value, const char *message)
{
    if (!value) { throw std::invalid_argument(message); }
}

void Finite(double value, const char *message)
{
    Require(std::isfinite(value), message);
}

double Cast(Real value)
{
    const double result = static_cast<double>(value);
    if (!std::isfinite(result))
    {
        throw std::range_error("Robust RANS result is not representable");
    }
    return result;
}

RANSClosureValues Modify(RANSClosureValues result, RANSModel model, double k,
    double epsilon, const RANSTensor &gradient, const RANSParameters &parameters,
    const RANSRobustVariant &variant)
{
    if (variant.Original()) { return result; }
    Require(model == RANSModel::CLSHighRe,
            "The declared robust constitutive variant requires CLSHighRe");
    Real production = 0;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            const Real increment = (static_cast<Real>(variant.beta)-1)*result.cubicStress[i][j]
                -static_cast<Real>(variant.gamma)*(static_cast<Real>(gradient[i][j])+
                    gradient[j][i]-variant.referenceStrain[i][j]);
            result.reynoldsStress[i][j] = Cast(static_cast<Real>(result.reynoldsStress[i][j])+increment);
            result.cubicStress[i][j] = Cast(static_cast<Real>(variant.beta)*result.cubicStress[i][j]);
            production -= static_cast<Real>(result.reynoldsStress[i][j])*gradient[i][j];
        }
    }
    result.production = Cast(production);
    result.kSource = Cast(production-epsilon);
    result.epsilonProduction = Cast(static_cast<Real>(parameters.cEpsilon1)*epsilon/k*production);
    result.epsilonSource = Cast(static_cast<Real>(result.epsilonProduction)-
        result.epsilonDestruction+result.yap);
    return result;
}

struct Limiter
{
    Real value, derivative;
};

Limiter Limit(Real s, Real ds, const RANSCappedDestructionControls &controls)
{
    const Real theta = controls.theta;
    if (s > theta || (s == theta && ds >= 0)) { return {1, 0}; }
    if (controls.smooth)
    {
        if (s <= 0) { return {0, 0}; }
        const Real t = s/theta, one = 1-t;
        return {t*t*t*(10+t*(-15+6*t)), 30*t*t*one*one*ds/theta};
    }
    return {s/theta, ds/theta};
}

void CheckPoint(double k, double K, double preE, double preD)
{
    Require(std::isfinite(k) && k >= 0, "Capped destruction requires finite nonnegative transported k");
    Require(std::isfinite(K) && K > 0, "Capped destruction requires positive finite sampled K");
    Require(std::isfinite(preE) && preE > 0, "Capped destruction requires positive finite pre-dissipation");
    Require(std::isfinite(preD) && preD >= 0, "Capped destruction requires nonnegative finite pre-destruction");
}

struct Budgets
{
    Real E = 0, D = 0, weighted = 0;
};

Budgets CheckInput(const RANSCappedDestructionInput &input,
                  const RANSCappedDestructionControls &controls)
{
    controls.Validate();
    const std::size_t n = input.k.size();
    Require(n > 0 && input.preE.size() == n && input.preD.size() == n &&
        input.metricWeights.size() == n, "Capped destruction input shapes differ or are empty");
    Budgets result;
    for (std::size_t i = 0; i < n; ++i)
    {
        CheckPoint(input.k[i], input.sampledK, input.preE[i], input.preD[i]);
        Require(std::isfinite(input.metricWeights[i]) && input.metricWeights[i] > 0,
                "Capped destruction metric weights must be positive and finite");
        const Real weight = input.metricWeights[i], destruction = weight*input.preD[i];
        result.E += weight*input.preE[i]; result.D += destruction;
        result.weighted += destruction*Limit(static_cast<Real>(input.k[i])/input.sampledK, 0, controls).value;
    }
    Require(result.E > 0 && result.D > 0, "Capped destruction requires a positive common pre-budget");
    Require(std::abs(result.E-result.D) <= controls.budgetRelativeTolerance*std::max(result.E,result.D),
            "Capped destruction pre-dissipation and pre-destruction budgets disagree");
    Cast(result.E); Cast(result.D);
    return result;
}

double Direction(const std::vector<double> &values, std::size_t i)
{
    return values.empty() ? 0 : values[i];
}

void CheckDirection(const RANSCappedDestructionDirection &direction, std::size_t size)
{
    Finite(direction.sampledK, "Nonfinite sampled-K direction");
    for (const auto *values : {&direction.preE,&direction.preD,&direction.k,&direction.metricWeights})
    {
        Require(values->empty() || values->size() == size, "Capped destruction direction shapes differ");
        for (double value : *values) { Finite(value, "Nonfinite capped-destruction direction"); }
    }
}

struct PointPair
{
    RANSCappedDestructionPoint value, derivative;
};

PointPair Point(double k, double K, double preE, double preD, Real b, Real m,
    Real dk, Real dK, Real dE, Real dD, Real db, Real dm,
    const RANSCappedDestructionControls &controls)
{
    const Real s = static_cast<Real>(k)/K;
    const Real ds = (dk-s*dK)/K;
    const auto limiter = Limit(s, ds, controls);
    const Real lambda = b/m, dlambda = (db-lambda*dm)/m;
    const Real E = lambda*preE, dEpos = dlambda*preE+lambda*dE;
    const Real D = preD*limiter.value/m;
    const Real dDpos = (dD*limiter.value+preD*limiter.derivative-D*dm)/m;
    const Real factor = controls.smooth ? 15.0L/8 : 1;
    const Real bound = factor*preD/(static_cast<Real>(controls.delta)*controls.theta*K);
    const Real dbound = factor*(dD-static_cast<Real>(preD)*dK/K)/
        (static_cast<Real>(controls.delta)*controls.theta*K);
    Real rate = 0, drate = 0;
    if (k > 0)
    {
        rate = D/k; drate = (dDpos-rate*dk)/k;
    }
    else if (!controls.smooth)
    {
        rate = static_cast<Real>(preD)/(m*controls.theta*K);
        drate = (dD-static_cast<Real>(preD)*(dm/m+dK/K))/(m*controls.theta*K);
    }
    return {{Cast(limiter.value),Cast(E),Cast(D),Cast(E-D),Cast(rate),Cast(bound)},
            {Cast(limiter.derivative),Cast(dEpos),Cast(dDpos),Cast(dEpos-dDpos),Cast(drate),Cast(dbound)}};
}
} // namespace

void RANSRobustVariant::Validate() const
{
    Require(std::isfinite(beta) && beta >= 0, "Common cubic beta must be finite and nonnegative");
    Require(std::isfinite(gamma) && gamma >= 0, "Prescribed gamma must be finite and nonnegative");
    Real trace = 0, scale = 1;
    for (int i = 0; i < 3; ++i)for (int j = 0; j < 3; ++j)
    {
        Finite(referenceStrain[i][j], "Nonfinite prescribed reference strain");
        scale = std::max(scale,std::abs(static_cast<Real>(referenceStrain[i][j])));
        if (i == j) { trace += referenceStrain[i][j]; }
    }
    if (gamma > 0)
    {
        Require(std::abs(trace) <= 64*std::numeric_limits<double>::epsilon()*scale,
                "Prescribed reference strain must be trace-free");
        for (int i = 0; i < 3; ++i)for (int j = 0; j < 3; ++j)
        {
            Require(std::abs(static_cast<Real>(referenceStrain[i][j])-referenceStrain[j][i]) <=
                64*std::numeric_limits<double>::epsilon()*scale,
                "Prescribed reference strain must be symmetric");
        }
    }
}

RANSClosureValues EvaluateRANSRobustClosure(RANSModel model, double k, double epsilon,
    const RANSTensor &gradient, double distance, const RANSParameters &parameters,
    double nu, const RANSRobustVariant &variant)
{
    variant.Validate();
    return Modify(EvaluateRANSClosure(model,k,epsilon,gradient,distance,parameters,nu),
                  model,k,epsilon,gradient,parameters,variant);
}

RANSClosureValues EvaluateRANSRobustClosureLinearization(RANSModel model, double k,
    double epsilon, const RANSTensor &gradient, double distance,
    const RANSParameters &parameters, double nu, bool strainBranch,
    const RANSRobustVariant &variant)
{
    variant.Validate();
    return Modify(EvaluateRANSClosureLinearization(model,k,epsilon,gradient,distance,parameters,nu,strainBranch),
                  model,k,epsilon,gradient,parameters,variant);
}

void RANSCappedDestructionControls::Validate() const
{
    Require(std::isfinite(theta) && theta > 0 && theta < 1, "Capped-destruction theta must lie in (0,1)");
    Require(std::isfinite(delta) && delta > 0 && delta < 1, "Capped-destruction delta must lie in (0,1)");
    Require(std::isfinite(budgetRelativeTolerance) && budgetRelativeTolerance > 0 &&
        budgetRelativeTolerance < 1, "Invalid capped-destruction common-budget tolerance");
}

RANSCappedDestructionPoint EvaluateRANSCappedDestructionPoint(double k, double K,
    double preE, double preD, double b, double m, const RANSCappedDestructionControls &controls)
{
    controls.Validate(); CheckPoint(k,K,preE,preD);
    Require(std::isfinite(b) && b >= 0 && b <= 1 && std::isfinite(m), "Invalid stored capped normalization");
    const double expected = std::max(b,controls.delta);
    Require(m == expected, "Stored capped denominator differs from max(b,delta)");
    return Point(k,K,preE,preD,b,m,0,0,0,0,0,0,controls).value;
}

RANSCappedDestructionValues EvaluateRANSCappedDestruction(const RANSCappedDestructionInput &input,
    const RANSCappedDestructionControls &controls)
{
    const auto budgets = CheckInput(input,controls);
    RANSCappedDestructionValues result;
    result.sampledK = input.sampledK; result.budgetE = Cast(budgets.E); result.budgetD = Cast(budgets.D);
    result.budgetMismatch = Cast(budgets.E-budgets.D);
    const Real b = budgets.weighted/budgets.D, m = std::max(b,static_cast<Real>(controls.delta));
    result.b = Cast(b); result.m = Cast(m); result.lambda = Cast(b/m); result.capped = b < controls.delta;
    Real E = 0, D = 0, C = 0;
    for (std::size_t i = 0; i < input.k.size(); ++i)
    {
        const auto value = Point(input.k[i],input.sampledK,input.preE[i],input.preD[i],b,m,0,0,0,0,0,0,controls).value;
        result.limiter.push_back(value.limiter); result.dissipation.push_back(value.dissipation);
        result.destruction.push_back(value.destruction); result.redistribution.push_back(value.redistribution);
        result.rate.push_back(value.rate); result.rateBound.push_back(value.rateBound);
        E += static_cast<Real>(input.metricWeights[i])*value.dissipation;
        D += static_cast<Real>(input.metricWeights[i])*value.destruction;
        C += static_cast<Real>(input.metricWeights[i])*value.redistribution;
    }
    result.dissipationBudget = Cast(E); result.destructionBudget = Cast(D); result.redistributionBudget = Cast(C);
    return result;
}

RANSCappedDestructionPoint EvaluateRANSCappedDestructionPoint(double preE, double preD,
    double localK, double K, const RANSCappedDestructionValues &normalization,
    const RANSCappedDestructionControls &controls)
{
    return EvaluateRANSCappedDestructionPoint(localK,K,preE,preD,
        normalization.b,normalization.m,controls);
}

RANSCappedDestructionValues RANSCappedDestructionDirectional(const RANSCappedDestructionInput &input,
    const RANSCappedDestructionDirection &direction, const RANSCappedDestructionControls &controls)
{
    const auto budgets = CheckInput(input,controls); CheckDirection(direction,input.k.size());
    Real dE = 0, dD = 0, dWeighted = 0;
    for (std::size_t i = 0; i < input.k.size(); ++i)
    {
        const Real dk = Direction(direction.k,i), dK = direction.sampledK;
        Require(input.k[i] != 0 || dk >= 0, "Capped-destruction direction leaves the nonnegative k cone");
        const Real weight = input.metricWeights[i], dw = Direction(direction.metricWeights,i);
        const Real D = input.preD[i], dlocalD = Direction(direction.preD,i);
        const Real s = static_cast<Real>(input.k[i])/input.sampledK;
        const auto limiter = Limit(s,(dk-s*dK)/input.sampledK,controls);
        dE += dw*input.preE[i]+weight*Direction(direction.preE,i);
        dD += dw*D+weight*dlocalD;
        dWeighted += (dw*D+weight*dlocalD)*limiter.value+weight*D*limiter.derivative;
    }
    const Real b = budgets.weighted/budgets.D, db = (dWeighted-b*dD)/budgets.D;
    const Real m = std::max(b,static_cast<Real>(controls.delta));
    const Real dm = b > controls.delta ? db : b < controls.delta ? 0 : std::max(db,0.0L);
    RANSCappedDestructionValues result;
    result.sampledK = direction.sampledK; result.budgetE = Cast(dE); result.budgetD = Cast(dD);
    result.budgetMismatch = Cast(dE-dD); result.b = Cast(db); result.m = Cast(dm);
    result.lambda = Cast((db-(b/m)*dm)/m); result.capped = b < controls.delta;
    Real E = 0, D = 0, C = 0;
    for (std::size_t i = 0; i < input.k.size(); ++i)
    {
        const auto value = Point(input.k[i],input.sampledK,input.preE[i],input.preD[i],b,m,
            Direction(direction.k,i),direction.sampledK,Direction(direction.preE,i),Direction(direction.preD,i),db,dm,controls);
        result.limiter.push_back(value.derivative.limiter); result.dissipation.push_back(value.derivative.dissipation);
        result.destruction.push_back(value.derivative.destruction); result.redistribution.push_back(value.derivative.redistribution);
        result.rate.push_back(value.derivative.rate); result.rateBound.push_back(value.derivative.rateBound);
        const Real weight = input.metricWeights[i], dw = Direction(direction.metricWeights,i);
        E += weight*value.derivative.dissipation+dw*value.value.dissipation;
        D += weight*value.derivative.destruction+dw*value.value.destruction;
        C += weight*value.derivative.redistribution+dw*value.value.redistribution;
    }
    result.dissipationBudget = Cast(E); result.destructionBudget = Cast(D); result.redistributionBudget = Cast(C);
    return result;
}
} // namespace Nektar::SolverUtils
