///////////////////////////////////////////////////////////////////////////////
// File: EquilibriumAWFProfiles.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <SolverUtils/Turbulence/EquilibriumAWFProfiles.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace Nektar::SolverUtils
{
namespace
{
void Finite(double x)
{
    if (!std::isfinite(x)) throw std::invalid_argument("Nonfinite EP-AWF input");
}
void Positive(double x)
{
    Finite(x);
    if (!(x > 0)) throw std::invalid_argument("EP-AWF requires positive input");
}
double Cast(long double x)
{
    const double result = static_cast<double>(x);
    if (!std::isfinite(result)) throw std::range_error("Unrepresentable EP-AWF result");
    return result;
}
AWFGeometry LegacyGeometry(const AWFProfileGeometry &g)
{
    g.Validate();
    // Original() needs a historical query point for its legacy value bundle.
    // That unused query is internal, fixed at H/2. It defines neither r_K nor
    // r_m; no legacy point dissipation or P/N interpolation enters this class.
    return {g.molecularNu, g.height, g.roughness, g.height/2};
}
double Gate(double x, double lo, double hi)
{
    if (x <= lo) return 1;
    if (x >= hi) return 0;
    const long double s = (static_cast<long double>(x)-lo)/(hi-lo), u = 1-s;
    return Cast(u*u*u*(1+3*s+6*s*s));
}
}

void AWFProfileGeometry::Validate() const
{
    Positive(molecularNu); Positive(height); Finite(roughness);
    if (!(roughness >= 0 && roughness < height && sampleHeight > 0 &&
          sampleHeight <= height && tailStart > 0 && tailStart <= height))
        throw std::invalid_argument("EP-AWF requires 0<=h<H and 0<r_K,r_m<=H");
    Finite(sampleHeight); Finite(tailStart);
}
double AWFProfileGeometry::SampleK(const std::function<double(double)> &field) const
{
    Validate();
    const double result = field(sampleHeight);
    Positive(result);
    return result;
}
void AWFNormalMetric::Validate(double height) const
{
    Positive(height); Finite(kappa1); Finite(kappa2); Positive(minimumStretch);
    if (!(minimumStretch <= 1 &&
          1+std::min(0.0, height*kappa1) >= minimumStretch &&
          1+std::min(0.0, height*kappa2) >= minimumStretch))
        throw std::invalid_argument("Invalid EP-AWF normal-ray metric");
    Positive(Jacobian(height));
}
double AWFNormalMetric::Jacobian(double r) const
{
    return Cast((1+static_cast<long double>(kappa1)*r)*
                (1+static_cast<long double>(kappa2)*r));
}
void AWFProfileSwitch::Validate() const
{
    Positive(on); Positive(off); Positive(marginOff); Positive(marginOn);
    Positive(reynoldsOff); Positive(reynoldsOn);
    if (!(on < off && marginOff < marginOn && reynoldsOff < reynoldsOn))
        throw std::invalid_argument("Unordered EP-AWF gate thresholds");
    for (double w : weights) Positive(w);
}
double EquilibriumAWFProfileWeight(double spatial, const std::array<double, 7> &sensors,
    double margin, double reynolds, const AWFProfileSwitch &c)
{
    c.Validate(); Finite(spatial); Finite(margin); Finite(reynolds);
    if (!(spatial >= 0 && spatial <= 1 && reynolds >= 0))
        throw std::invalid_argument("Invalid EP-AWF spatial/Reynolds weight");
    if (spatial == 0 || margin <= c.marginOff || reynolds <= c.reynoldsOff) return 0;
    long double d2 = 0;
    for (int j = 0; j < 7; ++j)
    {
        Finite(sensors[j]);
        d2 += static_cast<long double>(c.weights[j])*sensors[j]*sensors[j];
    }
    if (d2 >= static_cast<long double>(c.off)*c.off) return 0;
    const double state = Gate(Cast(d2), c.on*c.on, c.off*c.off);
    return spatial*(1-Gate(margin, c.marginOff, c.marginOn))*
                   (1-Gate(reynolds, c.reynoldsOff, c.reynoldsOn))*state;
}

AWFParameters CalibratedAWFParameters(const RANSNeutralReference &model)
{
    if (model.Model() != RANSModel::CLSHighRe)
        throw std::invalid_argument("CLS-calibrated AWF requires a CLS reference");
    AWFParameters wall;
    const long double q = model.Q();
    wall.cMu = Cast(1/(q*q));
    wall.cLength = Cast(static_cast<long double>(model.Kappa())*q*std::sqrt(q));
    wall.Validate();
    return wall;
}

AWFParameters WJCalibratedAWFParameters(const RANSNeutralReference &model)
{
    if (model.Model()!=RANSModel::WJHighRe)
        throw std::invalid_argument("WJ calibration requires a WJ reference");
    return {model.Cmu(), Cast(static_cast<long double>(model.Kappa())*
                std::pow(static_cast<long double>(model.Q()),1.5L)), 10.7};
}
EquilibriumAWFProfiles EquilibriumAWFProfiles::WJ(const AWFProfileGeometry &geometry,
    double k, const RANSNeutralReference &model, const AWFProfileOptions &options,
    const AWFNormalMetric &metric)
{
    if (options.formulation!=AWFProfileFormulation::WJCalibratedEPAWF &&
        options.formulation!=AWFProfileFormulation::WJCalibratedBudgetAWF)
        throw std::invalid_argument("WJ wall factory requires an explicit WJ source policy");
    return {geometry,k,model,metric,options,WJCalibratedAWFParameters(model)};
}

EquilibriumAWFProfiles::EquilibriumAWFProfiles(const AWFProfileGeometry &geometry, double k,
    const RANSNeutralReference &model, const AWFNormalMetric &metric,
    const AWFProfileOptions &options, const AWFParameters &wall)
    : m_geometry(geometry), m_metric(metric), m_options(options),
      m_sigmaK(model.Parameters().sigmaK),
      m_original(LegacyGeometry(geometry), k, wall)
{
    metric.Validate(geometry.height);
    Finite(options.requestedWeight); Positive(options.minimumReynolds);
    if (!(options.requestedWeight >= 0 && options.requestedWeight <= 1 &&
          options.quadratureOrder >= 2 && options.quadratureOrder <= 128))
        throw std::invalid_argument("Invalid EP-AWF weight or quadrature order");
    if (options.failurePolicy != AWFFailurePolicy::UseOriginal &&
        options.failurePolicy != AWFFailurePolicy::RequirePreservation)
        throw std::invalid_argument("Unknown EP-AWF failure policy");
    if (options.formulation != AWFProfileFormulation::LegacyEPAWF && !AcceptedAWF() && !WJCalibrated())
        throw std::invalid_argument("Unknown AWF profile formulation");
    if ((model.Model()==RANSModel::WJHighRe)!=WJCalibrated())
        throw std::invalid_argument("WJ outer reference requires an explicit WJ wall selection");
    if (AcceptedAWF() || WJCalibrated())
    {
        const auto calibrated = WJCalibrated() ? WJCalibratedAWFParameters(model) : CalibratedAWFParameters(model);
        const auto same = [](double a, double b) {
            return std::abs(a-b) <= 64*std::numeric_limits<double>::epsilon()*
                                      std::max(std::abs(a), std::abs(b));
        };
        if (geometry.sampleHeight != geometry.height || (AcceptedAWF() && geometry.tailStart != geometry.height))
            throw std::invalid_argument("CLSCalibratedBudgetAWF requires endpoint z_K=z_m=H");
        if (AcceptedAWF() && (options.requestedWeight != 1 ||
            options.failurePolicy != AWFFailurePolicy::RequirePreservation))
            throw std::invalid_argument("CLSCalibratedBudgetAWF has no chi or fallback policy");
        if (!same(wall.cMu, calibrated.cMu) || !same(wall.cLength, calibrated.cLength) ||
            wall.smoothIntercept != calibrated.smoothIntercept)
            throw std::invalid_argument("CLSCalibratedBudgetAWF requires explicit calibrated wall parameters");
    }
    for (double r : options.physicalBreakpoints)
    {
        Finite(r);
        if (!(r >= 0 && r <= geometry.height))
            throw std::invalid_argument("EP-AWF physical breakpoint outside wall layer");
    }
    GaussRule();
    auto &ref = m_reference;
    ref.requestedWeight = options.requestedWeight;
    ref.frictionVelocity = std::sqrt(k/model.Q());
    ref.traction = k/model.Q();
    ref.matchingVelocity = Cast(static_cast<long double>(ref.traction)*m_original.Integrals().j0);
    ref.numerator = Cast(std::pow(static_cast<long double>(ref.frictionVelocity), 3)/model.Kappa());
    Positive(ref.numerator);
    ref.reynolds = Cast(static_cast<long double>(geometry.height)*std::sqrt(k)/geometry.molecularNu);
    if (WJCalibrated() && !AcceptedAWF() && options.requestedWeight==0)
    {
        // Direct calibrated baseline: never evaluate E_c or its gate/root.
        m_quadrature=MakeQuadrature(geometry.height,0);
        return;
    }
    if (AcceptedAWF() || WJCalibrated())
    {
        // The new model has a budget applicability condition, not the legacy
        // chi/transition Reynolds gate. Re remains a reported scale only.
        // On a flat column the native closed AWF moment avoids cancellation
        // in J0-nu*J2 and fixes the counterfactual production budget directly.
        ref.productionBudget = metric.kappa1 == 0 && metric.kappa2 == 0 ?
            Cast(static_cast<long double>(ref.traction)*ref.traction*
                 m_original.Integrals().productionWeight) :
            Integrate([&](double r) { return ReferenceProduction(r); }, geometry.height, 0);
        ref.minimumBudget = CorrectedIntegral(geometry.height, geometry.tailStart);
        ref.margin = Cast((static_cast<long double>(ref.productionBudget)-ref.minimumBudget)/ref.minimumBudget);
        const bool lowRe = !AcceptedAWF() && ref.reynolds < options.minimumReynolds;
        if (ref.productionBudget < ref.minimumBudget || lowRe)
        {
            if (AcceptedAWF() || options.failurePolicy==AWFFailurePolicy::RequirePreservation)
                throw std::domain_error("Calibrated AWF cutoff/regime inadmissible; no cutoff clipping");
            ref.reason=lowRe ? AWFFallbackReason::OutsideHighReRegime : AWFFallbackReason::CutoffInadmissible;
            m_quadrature=MakeQuadrature(geometry.height,0);
            return; // Explicitly selected direct-baseline fallback, with weight=0.
        }
        if (metric.kappa1 == 0 && metric.kappa2 == 0)
        {
            const long double lambda = static_cast<long double>(ref.productionBudget)/ref.numerator;
            ref.cutoff = Cast(static_cast<long double>(geometry.height)*std::exp(1-lambda));
        }
        else
        {
            // The explicit metric-weighted dissipation integral is strictly
            // decreasing in c. A log bracket retains positive trial cutoffs.
            const double highStart = std::log(geometry.tailStart);
            double smallest = std::max(std::numeric_limits<double>::denorm_min(),
                                      ref.numerator/std::numeric_limits<double>::max());
            if (static_cast<long double>(ref.numerator)/smallest > std::numeric_limits<double>::max())
                smallest = std::nextafter(smallest, std::numeric_limits<double>::infinity());
            double lowLimit = std::log(smallest);
            while (std::exp(lowLimit) == 0 ||
                   static_cast<long double>(ref.numerator)/std::exp(lowLimit) > std::numeric_limits<double>::max())
                lowLimit = std::nextafter(lowLimit, std::numeric_limits<double>::infinity());
            double low = highStart, high = highStart, span = 1;
            while (BudgetDifference(std::exp(low)) < 0 && low > lowLimit)
            { low = std::max(lowLimit, highStart-span); span *= 2; }
            if (BudgetDifference(std::exp(low)) < 0)
                throw std::range_error("CLSCalibratedBudgetAWF cutoff_unrepresentable");
            for (int it = 0; it < 192; ++it)
            {
                const double mid = low+(high-low)/2;
                if (mid == low || mid == high) break;
                ++ref.cutoffIterations;
                if (BudgetDifference(std::exp(mid)) > 0) low = mid;
                else high = mid;
            }
            ref.cutoff = std::exp(low+(high-low)/2);
        }
        if (!(ref.cutoff > 0 && ref.cutoff <= geometry.tailStart))
            throw std::range_error("CLSCalibratedBudgetAWF cutoff_unrepresentable");
        Corrected(0, ref.cutoff); // Reject an unrepresentable plateau, never floor it.
        ref.weight = AcceptedAWF() ? 1 : options.requestedWeight; // Budget has no physical chi.
        ref.reason = AWFFallbackReason::Active;
        ref.budgetResidual = BudgetDifference(ref.cutoff);
        m_quadrature = MakeQuadrature(geometry.height, ref.cutoff);
        return;
    }
    ref.productionBudget = Integrate([&](double r) { return ReferenceProduction(r); },
                                     geometry.height, 0);
    ref.minimumBudget = Integrate([&](double r) { return Corrected(r, geometry.tailStart); },
                                  geometry.height, 0);
    Positive(ref.minimumBudget);
    ref.margin = Cast((static_cast<long double>(ref.productionBudget)-ref.minimumBudget)/ref.minimumBudget);
    auto fallback = [&](AWFFallbackReason reason) {
        ref.reason = reason;
        if (options.failurePolicy == AWFFailurePolicy::RequirePreservation)
            throw std::domain_error("Requested EP-AWF profile is inadmissible; no cutoff clipping");
    };
    if (options.requestedWeight > 0)
    {
        if (ref.reynolds < options.minimumReynolds)
            fallback(AWFFallbackReason::OutsideHighReRegime);
        else if (ref.productionBudget < ref.minimumBudget)
            fallback(AWFFallbackReason::CutoffInadmissible);
        else
        {
          try
          {
            // A bracket in log(r_c) preserves strictly positive trial cutoffs.
            const double upper = std::log(geometry.tailStart);
            // Subnormal cutoffs can be valid. Bound the search by a
            // representable positive cutoff AND finite A/r_c, rather than
            // the smallest normal number. Correct log/exp rounding upward
            // at this numerical range boundary before evaluating a trial.
            const long double maximum = std::numeric_limits<double>::max();
            double smallest = Cast(std::max(
                static_cast<long double>(std::numeric_limits<double>::denorm_min()),
                static_cast<long double>(ref.numerator)/maximum));
            if (static_cast<long double>(ref.numerator)/smallest > maximum)
                smallest = std::nextafter(smallest, std::numeric_limits<double>::infinity());
            double lowerLimit = std::log(smallest);
            while (std::exp(lowerLimit) == 0 ||
                   static_cast<long double>(ref.numerator)/std::exp(lowerLimit) > maximum)
                lowerLimit = std::nextafter(lowerLimit, std::numeric_limits<double>::infinity());
            double low = upper, high = upper, span = 1;
            double residual = BudgetDifference(std::exp(low));
            while (residual < 0 && low > lowerLimit)
            {
                low = std::max(lowerLimit, upper-span);
                residual = BudgetDifference(std::exp(low));
                span *= 2;
            }
            if (residual < 0)
                fallback(AWFFallbackReason::UnrepresentableCutoff);
            else
            {
                for (int it = 0; it < 192; ++it)
                {
                    const double mid = low+(high-low)/2;
                    if (mid == low || mid == high) break;
                    ++ref.cutoffIterations;
                    if (BudgetDifference(std::exp(mid)) > 0) low = mid;
                    else high = mid;
                }
                ref.cutoff = std::exp(low+(high-low)/2);
                if (!(ref.cutoff > 0 && ref.cutoff <= geometry.tailStart))
                    fallback(AWFFallbackReason::UnrepresentableCutoff);
                else
                {
                    // Also reject an unrepresentable cap; never hide it with
                    // a dissipation floor or a clipped cutoff.
                    Corrected(0, ref.cutoff);
                    ref.weight = options.requestedWeight;
                    ref.reason = AWFFallbackReason::Active;
                }
            }
          }
          catch (const std::range_error &)
          {
              ref.cutoff = 0;
              ref.weight = 0;
              fallback(AWFFallbackReason::UnrepresentableCutoff);
          }
        }
    }
    m_quadrature = MakeQuadrature(geometry.height, ref.weight > 0 ? ref.cutoff : 0);
    // Form both discrete budgets on the final source quadrature.
    ref.productionBudget = Integrate([&](double r) { return ReferenceProduction(r); },
                                     geometry.height, ref.weight > 0 ? ref.cutoff : 0);
    ref.minimumBudget = Integrate([&](double r) { return Corrected(r, geometry.tailStart); },
                                  geometry.height, ref.weight > 0 ? ref.cutoff : 0);
    ref.margin = Cast((static_cast<long double>(ref.productionBudget)-ref.minimumBudget)/ref.minimumBudget);
    if (ref.weight > 0) ref.budgetResidual = BudgetDifference(ref.cutoff);
}

void EquilibriumAWFProfiles::CheckHeight(double r) const
{
    Finite(r);
    if (!(r >= 0 && r <= m_geometry.height))
        throw std::invalid_argument("EP-AWF evaluation outside [0,H]");
}
EquilibriumAWFProfiles EquilibriumAWFProfiles::Reweighted(double weight) const
{
    if (AcceptedAWF())
        throw std::invalid_argument("CLSCalibratedBudgetAWF has no chi reweighting");
    Finite(weight);
    if (!(weight>0 && weight<=1 && m_reference.weight>0))
        throw std::invalid_argument("Reweighting requires an admissible active cutoff and 0<chi<=1");
    auto result=*this;
    result.m_options.requestedWeight=weight;
    result.m_reference.requestedWeight=weight;
    result.m_reference.weight=weight;
    return result;
}
void EquilibriumAWFProfiles::GaussRule()
{
    const int n = m_options.quadratureOrder;
    m_gaussPoints.resize(n); m_gaussWeights.resize(n);
    for (int j = 0; j < (n+1)/2; ++j)
    {
        long double x = std::cos(std::acos(-1.0L)*(j+.75L)/(n+.5L)), derivative = 0;
        for (int it = 0; it < 100; ++it)
        {
            long double p = 1, previous = 0;
            for (int k = 1; k <= n; ++k)
            {
                const long double next = ((2*k-1)*x*p-(k-1)*previous)/k;
                previous = p; p = next;
            }
            derivative = n*(x*p-previous)/(x*x-1);
            const long double delta = p/derivative;
            x -= delta;
            if (std::abs(delta) <= 4*std::numeric_limits<long double>::epsilon()) break;
        }
        // Re-evaluate at the final root before forming its positive weight.
        long double p = 1, previous = 0;
        for (int k = 1; k <= n; ++k)
        {
            const long double next = ((2*k-1)*x*p-(k-1)*previous)/k;
            previous = p; p = next;
        }
        derivative = n*(x*p-previous)/(x*x-1);
        const double weight = Cast(2/((1-x*x)*derivative*derivative));
        Positive(weight);
        m_gaussPoints[j] = Cast(-x); m_gaussPoints[n-1-j] = Cast(x);
        m_gaussWeights[j] = weight; m_gaussWeights[n-1-j] = weight;
    }
}
std::vector<AWFProfileQuadraturePoint> EquilibriumAWFProfiles::MakeQuadrature(
    double upper, double cutoff, double lower) const
{
    CheckHeight(upper); CheckHeight(lower);
    if (lower > upper) throw std::invalid_argument("Reversed AWF integration interval");
    if (upper == lower) return {};
    std::vector<double> breaks{lower, upper};
    for (double r : m_options.physicalBreakpoints)
        if (r > lower && r < upper) breaks.push_back(r);
    for (double r : {m_geometry.roughness, m_original.ViscosityInterceptHeight(),
                     m_original.DissipationCutoff(), m_geometry.tailStart, cutoff})
        if (r > lower && r < upper) breaks.push_back(r);
    const auto &p = m_original.Parameters();
    long double r = std::min(static_cast<long double>(m_geometry.tailStart),
        static_cast<long double>(m_geometry.molecularNu)/
        (p.cMu*p.cLength*std::sqrt(static_cast<long double>(K()))));
    if (cutoff > 0) r = std::min(r, static_cast<long double>(cutoff));
    for (int j = 0; j < 4096 && r < upper; ++j, r *= 2)
        if (r > lower && Cast(r) > lower) breaks.push_back(Cast(r));
    std::sort(breaks.begin(), breaks.end());
    breaks.erase(std::unique(breaks.begin(), breaks.end()), breaks.end());
    std::vector<AWFProfileQuadraturePoint> result;
    for (std::size_t j = 1; j < breaks.size(); ++j)
    {
        const long double half = (static_cast<long double>(breaks[j])-breaks[j-1])/2;
        const long double mid = breaks[j-1]+half;
        for (int k = 0; k < m_options.quadratureOrder; ++k)
        {
            const double point = Cast(mid+half*m_gaussPoints[k]);
            const double weight = Cast(half*m_gaussWeights[k]);
            Positive(weight);
            const double jacobian = m_metric.Jacobian(point); Positive(jacobian);
            result.push_back({point, weight, jacobian});
        }
    }
    return result;
}
std::vector<AWFProfileQuadraturePoint> EquilibriumAWFProfiles::IntegrationQuadrature(
    double upper,const std::vector<double> &additional) const
{
    auto copy=*this;
    for(double z:additional){Finite(z);if(z>0 && z<m_geometry.height)copy.m_options.physicalBreakpoints.push_back(z);}
    return copy.MakeQuadrature(upper,m_reference.weight>0?m_reference.cutoff:0.);
}
double EquilibriumAWFProfiles::Integrate(const std::function<double(double)> &f,
    double upper, double cutoff) const
{
    long double result = 0;
    for (const auto &q : MakeQuadrature(upper, cutoff))
        result += static_cast<long double>(q.weight)*q.jacobian*f(q.r);
    return Cast(result);
}
double EquilibriumAWFProfiles::Corrected(double r, double cutoff) const
{
    CheckHeight(r); Positive(cutoff);
    return Cast(static_cast<long double>(m_reference.numerator)/std::max(r, cutoff));
}
double EquilibriumAWFProfiles::CorrectedIntegral(double upper, double cutoff) const
{
    CheckHeight(upper); Positive(cutoff);
    const long double r = upper, c = cutoff;
    const long double s = static_cast<long double>(m_metric.kappa1)+m_metric.kappa2;
    const long double p = static_cast<long double>(m_metric.kappa1)*m_metric.kappa2;
    const long double integral = r <= c ? (r+s*r*r/2+p*r*r*r/3)/c :
        1+std::log(r/c)+s*(r-c/2)+p*(r*r/2-c*c/6);
    return Cast(static_cast<long double>(m_reference.numerator)*integral);
}
double EquilibriumAWFProfiles::DissipationBudget(double cutoff) const
{
    Positive(cutoff);
    if (cutoff > m_geometry.height)
        throw std::invalid_argument("AWF budget cutoff outside (0,H]");
    return CorrectedIntegral(m_geometry.height, cutoff);
}
double EquilibriumAWFProfiles::BudgetDifference(double cutoff) const
{
    if (AcceptedAWF() || WJCalibrated())
        return Cast(static_cast<long double>(CorrectedIntegral(m_geometry.height, cutoff))-
                    m_reference.productionBudget);
    return Integrate([&](double r) { return Corrected(r, cutoff)-ReferenceProduction(r); },
                     m_geometry.height, cutoff);
}
double EquilibriumAWFProfiles::Production(double r, const AWFVector &u, const AWFVector &f) const
{
    CheckHeight(r);
    const auto derivative = m_original.VelocityGradient(r, u, f);
    return Cast(static_cast<long double>(m_original.EddyViscosity(r))*
                (static_cast<long double>(derivative[0])*derivative[0]+
                 static_cast<long double>(derivative[1])*derivative[1]));
}
double EquilibriumAWFProfiles::ReferenceProduction(double r) const
{
    if (AcceptedAWF())
    {
        CheckHeight(r);
        const long double d = m_original.TotalViscosity(r);
        return Cast(static_cast<long double>(m_reference.traction)*m_reference.traction*
                    m_original.EddyViscosity(r)/(d*d));
    }
    // Same arithmetic path as actual production at the reference, so its
    // all-mode cancellation is pointwise before any source quadrature.
    return Production(r, {m_reference.matchingVelocity, 0}, {0, 0});
}
AWFProfileValues EquilibriumAWFProfiles::AtProduction(double r, double production) const
{
    CheckHeight(r); Finite(production);
    if (AcceptedAWF())
    {
        const double e = Corrected(r, m_reference.cutoff);
        return {production, ReferenceProduction(r), e, 0,
                Cast(static_cast<long double>(production)-e)};
    }
    const double original = m_original.OriginalDissipation(r), p0 = ReferenceProduction(r);
    const double chi = m_reference.weight;
    if (chi == 0) return {production, p0, original, 0, Cast(static_cast<long double>(production)-original)};
    const double corrected = Corrected(r, m_reference.cutoff);
    return {production, p0, Cast((1-static_cast<long double>(chi))*original+chi*corrected),
        Cast(static_cast<long double>(chi)*(corrected-p0)),
        Cast(static_cast<long double>(production)-(1-static_cast<long double>(chi))*original-chi*p0)};
}
AWFProfileValues EquilibriumAWFProfiles::Evaluate(double r, const AWFVector &u, const AWFVector &f) const
{
    return AtProduction(r, Production(r, u, f));
}
double EquilibriumAWFProfiles::Dissipation(double r) const
{
    return AtProduction(r, 0).dissipation;
}
double EquilibriumAWFProfiles::RedistributionFlux(double r) const
{
    CheckHeight(r);
    if (AcceptedAWF()) return 0;
    if (m_reference.weight == 0 || r == 0) return 0;
    return Cast(static_cast<long double>(m_reference.weight)/m_metric.Jacobian(r)*
        Integrate([&](double z) { return Corrected(z, m_reference.cutoff)-ReferenceProduction(z); },
                  r, m_reference.cutoff));
}
double EquilibriumAWFProfiles::ReferenceKDiffusivity(double r) const
{
    CheckHeight(r);
    return Cast(static_cast<long double>(m_geometry.molecularNu)+m_original.EddyViscosity(r)/m_sigmaK);
}
double EquilibriumAWFProfiles::ReferenceKFlux(double r) const
{
    CheckHeight(r);
    if (!AcceptedAWF() || r == 0) return 0;
    const double production = m_metric.kappa1 == 0 && m_metric.kappa2 == 0 ?
        Cast(static_cast<long double>(m_reference.traction)*m_reference.traction*
             m_original.Integrals(r).productionWeight) :
        Integrate([&](double z) { return ReferenceProduction(z); }, r, m_reference.cutoff);
    return Cast(static_cast<long double>(CorrectedIntegral(r, m_reference.cutoff))-production);
}
double EquilibriumAWFProfiles::ReferenceKGradient(double r) const
{
    CheckHeight(r);
    if (!AcceptedAWF()) return 0;
    return Cast(static_cast<long double>(ReferenceKFlux(r))/
                (static_cast<long double>(m_metric.Jacobian(r))*ReferenceKDiffusivity(r)));
}
double EquilibriumAWFProfiles::ReferenceK(double r) const
{
    CheckHeight(r);
    if (!AcceptedAWF() || r == m_geometry.height) return K();
    long double integral = 0;
    for (const auto &q : MakeQuadrature(m_geometry.height, m_reference.cutoff, r))
        integral += static_cast<long double>(q.weight)*ReferenceKGradient(q.r);
    return Cast(static_cast<long double>(K())-integral);
}
std::vector<double> EquilibriumAWFProfiles::RayLoads(
    const std::vector<std::vector<double>> &basis, const std::vector<double> &samples) const
{
    if (samples.size() != m_quadrature.size())
        throw std::invalid_argument("EP-AWF source samples do not match split quadrature");
    std::vector<double> result;
    for (const auto &row : basis)
    {
        if (row.size() != samples.size()) throw std::invalid_argument("Invalid EP-AWF basis row");
        long double load = 0;
        for (std::size_t j = 0; j < row.size(); ++j)
        {
            Finite(row[j]); Finite(samples[j]);
            load += static_cast<long double>(m_quadrature[j].weight)*m_quadrature[j].jacobian*
                    row[j]*samples[j];
        }
        result.push_back(Cast(load));
    }
    return result;
}
std::vector<double> EquilibriumAWFProfiles::SourceLoads(
    const std::vector<std::vector<double>> &basis, const AWFVector &u, const AWFVector &f) const
{
    std::vector<double> source;
    for (const auto &q : m_quadrature) source.push_back(Evaluate(q.r, u, f).kSource);
    return RayLoads(basis, source);
}
}
