///////////////////////////////////////////////////////////////////////////////
// File: AnalyticalWallFunction.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <SolverUtils/Turbulence/AnalyticalWallFunction.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace Nektar::SolverUtils
{
namespace
{
using Real = long double;
void Finite(double value,const char *name)
{
    if (!std::isfinite(value)) { throw std::invalid_argument(name); }
}
void Positive(double value,const char *name)
{
    Finite(value,name);
    if (!(value>0.0)) { throw std::invalid_argument(name); }
}
double Cast(Real value)
{
    const double result=static_cast<double>(value);
    if (!std::isfinite(result)) { throw std::range_error("AWF result is not representable"); }
    return result;
}
void Check(const AWFVector &vector)
{
    for (double v : vector) { Finite(v,"Nonfinite AWF tangential input"); }
}
Real LogRatio(Real x) { return x==0.0L ? 1.0L : std::log1p(x)/x; }
Real FirstMoment(Real x)
{
    if (x>0.5L) { return (1.0L-LogRatio(x))/x; }
    Real result=0.0L,power=1.0L;
    for (int n=0; n<160; ++n)
    {
        const Real term=power/(n+2.0L); result+=term;
        if (std::abs(term)<std::numeric_limits<Real>::epsilon()*std::abs(result)) { break; }
        power *= -x;
    }
    return result;
}
// A_j(x)=integral_0^1 t^j/(1+xt)^2 dt. The series avoids catastrophic
// subtraction near the molecular limit. The other expressions are exact.
Real SquaredMoment(int j,Real x)
{
    if (x<=0.5L)
    {
        Real result=0.0L,power=1.0L;
        for (int n=0; n<160; ++n)
        {
            const Real term=(n+1.0L)*power/(j+n+1.0L); result+=term;
            if (std::abs(term)<std::numeric_limits<Real>::epsilon()*std::abs(result)) { break; }
            power *= -x;
        }
        return result;
    }
    const Real inverse=1.0L/x, log=std::log1p(x), fraction=x/(1.0L+x);
    if (j==0) { return 1.0L/(1.0L+x); }
    if (j==1) { return (log-fraction)*inverse*inverse; }
    if (j==2) { return (1.0L+(-2.0L*log+fraction)*inverse)*inverse*inverse; }
    return (0.5L-2.0L*inverse+(3.0L*log-fraction)*inverse*inverse)*inverse*inverse;
}
} // namespace

void AWFGeometry::Validate() const
{
    Positive(molecularNu,"AWF requires positive physical molecular viscosity");
    Positive(height,"AWF requires positive layer height");
    Positive(pointP,"AWF requires positive P height");
    Finite(roughness,"Nonfinite roughness");
    if (!(roughness>=0.0 && roughness<height && pointP<height))
        { throw std::invalid_argument("AWF requires 0 <= h < H and 0 < P < H"); }
}
void AWFParameters::Validate() const
{
    for (double value : {cMu,cLength,smoothIntercept})
        { Positive(value,"AWF coefficients must be finite and positive"); }
}

AnalyticalWallFunction::AnalyticalWallFunction(
    const AWFGeometry &geometry,double k,const AWFParameters &p)
    : m_geometry(geometry),m_parameters(p),m_k(k)
{
    geometry.Validate(); p.Validate(); Positive(k,"AWF K must be finite and positive");
    const Real rootK=std::sqrt(static_cast<Real>(k)),nu=geometry.molecularNu;
    m_intercept=p.smoothIntercept;
    if (geometry.roughness>0.0)
    {
        const Real x=geometry.roughness*rootK/(70.0L*nu);
        const Real exponent=std::max(0.5L-0.4L*std::pow(x,0.7L),
                                    1.0L-0.79L*std::pow(x,-0.28L));
        m_intercept=-p.smoothIntercept*std::expm1(exponent*std::log(x));
    }
    m_zv=nu*m_intercept/rootK;
    m_slope=static_cast<Real>(p.cMu)*p.cLength*rootK;
    m_ze=2.0L*p.cLength*nu/rootK;
    Cast(m_intercept); Cast(m_zv); Positive(Cast(m_ze),"Unrepresentable AWF dissipation cutoff");
    m_integrals=Integrals(geometry.height);
    Positive(m_integrals.j0,"Unrepresentable AWF J0");
    Positive(MeanOriginalDissipation(),"Unrepresentable mean AWF dissipation");
}

void AnalyticalWallFunction::CheckHeight(double z) const
{
    Finite(z,"Nonfinite auxiliary height");
    if (!(z>=0.0 && z<=m_geometry.height)) { throw std::invalid_argument("Auxiliary height outside [0,H]"); }
}
double AnalyticalWallFunction::Intercept() const { return Cast(m_intercept); }
double AnalyticalWallFunction::ViscosityInterceptHeight() const { return Cast(m_zv); }
double AnalyticalWallFunction::DissipationCutoff() const { return Cast(m_ze); }
double AnalyticalWallFunction::EddyViscosity(double z) const
{
    CheckHeight(z); return Cast(m_slope*std::max(0.0L,static_cast<Real>(z)-m_zv));
}
double AnalyticalWallFunction::TotalViscosity(double z) const
{
    CheckHeight(z); return Cast(m_geometry.molecularNu+
                               m_slope*std::max(0.0L,static_cast<Real>(z)-m_zv));
}

std::array<Real,3> AnalyticalWallFunction::ProductionMoments(Real left,Real right) const
{
    const Real length=right-left;
    if (!(length>0.0L)) { return {}; }
    const Real nut=m_slope*std::max(0.0L,left-m_zv),D=m_geometry.molecularNu+nut;
    const Real x=m_slope*length/D;
    std::array<Real,3> moments{};
    Real power=length/(D*D);
    for (int j=0; j<3; ++j)
    {
        moments[j]=power*(nut*SquaredMoment(j,x)+m_slope*length*SquaredMoment(j+1,x));
        power*=length;
    }
    return moments;
}

AWFIntegrals AnalyticalWallFunction::Integrals(double z) const
{
    CheckHeight(z);
    const Real upper=z,nu=m_geometry.molecularNu,h=m_geometry.roughness;
    const Real split=std::min(upper,std::max(0.0L,m_zv));
    Real j0=split/nu,j2=split/(nu*nu),j1=std::pow(std::max(0.0L,split-h),2)/(2.0L*nu);
    Real weight=0.0L;
    if (upper>split)
    {
        const Real length=upper-split,D=nu+m_slope*(split-m_zv),x=m_slope*length/D;
        j0+=length/D*LogRatio(x);
        j2+=length/D/(D+m_slope*length);
        const Real start=std::max(split,h),span=upper-start;
        if (span>0.0L)
        {
            const Real Dl=nu+m_slope*(start-m_zv),ratio=m_slope*span/Dl;
            j1+=(start-h)*span/Dl*LogRatio(ratio)+span*span/Dl*FirstMoment(ratio);
        }
        weight=ProductionMoments(split,upper)[0];
    }
    return {Cast(j0),Cast(j1),Cast(j2),Cast(weight)};
}

double AnalyticalWallFunction::OriginalDissipation(double z) const
{
    CheckHeight(z);
    const Real numerator=static_cast<Real>(m_k)*std::sqrt(static_cast<Real>(m_k))/m_parameters.cLength;
    return Cast(numerator/std::max(static_cast<Real>(z),m_ze));
}
double AnalyticalWallFunction::MeanOriginalDissipation() const
{
    const Real numerator=static_cast<Real>(m_k)*std::sqrt(static_cast<Real>(m_k))/m_parameters.cLength;
    if (m_geometry.height<=m_ze) { return Cast(numerator/m_ze); }
    return Cast(numerator/m_geometry.height*(1.0L+std::log(m_geometry.height/m_ze)));
}

std::array<Real,2> AnalyticalWallFunction::Traction(const AWFVector &u,const AWFVector &f) const
{
    Check(u); Check(f);
    const Real B=static_cast<Real>(m_integrals.j1)/m_integrals.j0;
    return {u[0]/static_cast<Real>(m_integrals.j0)-B*f[0],
            u[1]/static_cast<Real>(m_integrals.j0)-B*f[1]};
}

AWFValues AnalyticalWallFunction::Evaluate(const AWFVector &u,const AWFVector &f) const
{
    const auto traction=Traction(u,f);
    const Real start=std::min(static_cast<Real>(m_geometry.height),std::max(0.0L,m_zv));
    const Real split=std::max(start,static_cast<Real>(m_geometry.roughness));
    const std::array<Real,3> boundaries{start,std::min(split,static_cast<Real>(m_geometry.height)),m_geometry.height};
    Real production=0.0L;
    for (int segment=0; segment<2; ++segment)
    {
        const Real left=boundaries[segment],right=boundaries[segment+1];
        if (!(right>left)) { continue; }
        const auto q=ProductionMoments(left,right);
        if (!(q[0]>0.0L)) { continue; }
        const Real mean=q[1]/q[0];
        Real variance=q[2]-q[1]*mean;
        const Real tolerance=64.0L*std::numeric_limits<Real>::epsilon()*std::max(q[2],q[1]*mean);
        if (variance < -tolerance) { throw std::range_error("Loss of positive AWF production moment"); }
        variance=std::max(0.0L,variance); // at most roundoff in a positive Gram matrix
        for (int component=0; component<2; ++component)
        {
            const Real g=left>=m_geometry.roughness ? f[component] : 0.0L;
            const Real T=traction[component]+f[component]*std::max(0.0L,left-m_geometry.roughness);
            const Real center=T+g*mean;
            production+=q[0]*center*center+variance*g*g;
        }
    }
    return {{Cast(traction[0]),Cast(traction[1])},Cast(production/m_geometry.height),
            MeanOriginalDissipation(),OriginalDissipation(m_geometry.pointP),
            OriginalDissipation(m_geometry.height)};
}

AWFVector AnalyticalWallFunction::Velocity(double z,const AWFVector &u,const AWFVector &f) const
{
    const auto traction=Traction(u,f); const auto integrals=Integrals(z);
    return {Cast(traction[0]*integrals.j0+static_cast<Real>(f[0])*integrals.j1),
            Cast(traction[1]*integrals.j0+static_cast<Real>(f[1])*integrals.j1)};
}
AWFVector AnalyticalWallFunction::VelocityGradient(double z,const AWFVector &u,const AWFVector &f) const
{
    CheckHeight(z); const auto traction=Traction(u,f);
    const Real height=std::max(0.0L,static_cast<Real>(z)-m_geometry.roughness);
    const Real D=m_geometry.molecularNu+m_slope*std::max(0.0L,static_cast<Real>(z)-m_zv);
    return {Cast((traction[0]+f[0]*height)/D),Cast((traction[1]+f[1]*height)/D)};
}

EquilibriumWallFunction::EquilibriumWallFunction(
    const AWFGeometry &geometry,double k,const RANSNeutralReference &model,const AWFParameters &p)
    : m_original(geometry,k,p),m_model(model),m_reference{}
{
    if (model.Model()==RANSModel::WJHighRe)
        throw std::invalid_argument("WJ requires the calibrated EquilibriumAWFProfiles::WJ factory");
    auto &r=m_reference;
    r.frictionVelocity=std::sqrt(k/model.Q()); r.traction=k/model.Q();
    r.matchingVelocity=Cast(static_cast<Real>(r.traction)*m_original.Integrals().j0);
    r.numerator=Cast(std::pow(static_cast<Real>(r.frictionVelocity),3)/model.Kappa());
    Positive(r.numerator,"Unrepresentable reference dissipation numerator");
    r.production=Cast(static_cast<Real>(r.traction)*r.traction/geometry.height*
                      m_original.Integrals().productionWeight);
    r.lambda=Cast(static_cast<Real>(geometry.height)*r.production/r.numerator);
    r.minimumLambda=Cast(1.0L+std::log(static_cast<Real>(geometry.height)/geometry.pointP));
    const Real logCutoff=std::log(static_cast<Real>(geometry.height))+1.0L-r.lambda;
    const Real cutoff=std::exp(logCutoff);
    r.cutoff=Cast(cutoff);
    r.logAerodynamicRoughness=Cast(std::log(static_cast<Real>(geometry.height))-
                                  model.Kappa()*r.matchingVelocity/r.frictionVelocity);
    r.reason=AWFFallbackReason::CutoffInadmissible;
    if (!(r.cutoff>0.0) || !std::isfinite(static_cast<double>(r.numerator/cutoff)))
        { r.reason=AWFFallbackReason::UnrepresentableCutoff; }
    else if (r.cutoff<=geometry.pointP)
        { r.admissible=true; r.reason=AWFFallbackReason::Active; }
}

void AWFSwitchParameters::Validate() const
{
    Positive(on,"Equilibrium switch requires on > 0"); Positive(off,"Invalid switch off threshold");
    if (!(off>on)) { throw std::invalid_argument("Equilibrium switch requires off > on"); }
    bool active=false;
    for (double weight : weights)
    {
        Finite(weight,"Nonfinite equilibrium sensor weight");
        if (weight<0.0) { throw std::invalid_argument("Negative equilibrium sensor weight"); }
        active|=weight>0.0;
    }
    Check(incomingDirection);
    if (!active || !(std::hypot(incomingDirection[0],incomingDirection[1])>0.0))
        { throw std::invalid_argument("Equilibrium switch requires active sensors and a fixed incoming direction"); }
}

double AWFEquilibriumWeight(double departure,double on,double off)
{
    Finite(departure,"Nonfinite equilibrium departure"); Positive(on,"Invalid on threshold"); Positive(off,"Invalid off threshold");
    if (departure<0.0 || !(off>on)) { throw std::invalid_argument("Invalid equilibrium switch arguments"); }
    if (departure<=on) { return 1.0; }
    if (departure>=off) { return 0.0; }
    const Real d=departure,a=on,b=off,s=((d-a)/(b-a))*((d+a)/(b+a)),u=1.0L-s;
    return Cast(u*u*u*(1.0L+3.0L*s+6.0L*s*s));
}

EquilibriumAWFValues EquilibriumWallFunction::EvaluateFixedWeight(
    const AWFVector &u,const AWFVector &f,double weight,AWFFailurePolicy policy,bool highRe) const
{
    Finite(weight,"Nonfinite equilibrium weight");
    if (weight<0.0 || weight>1.0) { throw std::invalid_argument("Equilibrium weight must lie in [0,1]"); }
    if (policy!=AWFFailurePolicy::UseOriginal && policy!=AWFFailurePolicy::RequirePreservation)
        { throw std::invalid_argument("Unknown AWF failure policy"); }
    EquilibriumAWFValues result{};
    result.values=m_original.Evaluate(u,f); result.requestedWeight=weight;
    result.reason=AWFFallbackReason::SwitchInactive;
    if (weight==0.0) { return result; } // exact complete baseline recovery
    result.reason=!highRe ? AWFFallbackReason::OutsideHighReRegime : m_reference.reason;
    if (!highRe || !m_reference.admissible)
    {
        if (policy==AWFFailurePolicy::RequirePreservation)
            { throw std::domain_error("Requested equilibrium AWF reference is not admissible in this wall regime"); }
        return result;
    }
    result.weight=weight; result.reason=AWFFallbackReason::Active;
    if (weight==1.0)
    {
        result.values.meanDissipation=m_reference.production;
    }
    else
    {
        result.values.meanDissipation=Cast((1.0L-weight)*result.values.meanDissipation+weight*static_cast<Real>(m_reference.production));
    }
    result.values.dissipationP=Dissipation(m_original.Geometry().pointP,weight);
    result.values.dissipationH=Dissipation(m_original.Geometry().height,weight);
    return result;
}

double EquilibriumWallFunction::Dissipation(double z,double weight) const
{
    const double original=m_original.OriginalDissipation(z);
    Finite(weight,"Nonfinite dissipation blend weight");
    if (weight<0.0 || weight>1.0) { throw std::invalid_argument("Invalid dissipation blend weight"); }
    if (weight==0.0) { return original; }
    if (!m_reference.admissible) { throw std::domain_error("Corrected dissipation is not admissible"); }
    const Real corrected=m_reference.numerator/std::max(static_cast<Real>(z),static_cast<Real>(m_reference.cutoff));
    return Cast((1.0L-weight)*original+weight*corrected);
}

EquilibriumAWFValues EquilibriumWallFunction::Evaluate(
    const AWFIndependentSamples &sample,const AWFSwitchParameters &p,AWFFailurePolicy policy,bool highRe) const
{
    p.Validate(); Check(sample.matchingVelocity); Check(sample.velocityP); Check(sample.velocityN);
    Check(sample.pressureForcing); Check(sample.convectionForcing);
    const auto &g=m_original.Geometry();
    Positive(sample.heightN,"Invalid independent N height");
    if (!(sample.heightN>g.height)) { throw std::invalid_argument("N must be above the matching face"); }
    for (double value : {sample.kN,sample.epsilonN,sample.normalVelocityP}) { Finite(value,"Nonfinite independent wall sample"); }
    const Real norm=std::hypot(static_cast<Real>(p.incomingDirection[0]),static_cast<Real>(p.incomingDirection[1]));
    const Real us=std::hypot(static_cast<Real>(m_reference.frictionVelocity),static_cast<Real>(g.molecularNu)/g.height);
    Real du=0.0L,delta=0.0L,fp=0.0L,fa=0.0L;
    AWFVector forcing;
    for (int i=0; i<2; ++i)
    {
        const Real direction=p.incomingDirection[i]/norm;
        du=std::hypot(du,(sample.matchingVelocity[i]-m_reference.matchingVelocity*direction)/us);
        delta=std::hypot(delta,(sample.velocityN[i]-sample.velocityP[i]-
            m_reference.frictionVelocity/m_model.Kappa()*std::log(static_cast<Real>(sample.heightN)/g.pointP)*direction)/us);
        fp=std::hypot(fp,g.height*static_cast<Real>(sample.pressureForcing[i])/(us*us));
        fa=std::hypot(fa,g.height*static_cast<Real>(sample.convectionForcing[i])/(us*us));
        forcing[i]=Cast(static_cast<Real>(sample.pressureForcing[i])+sample.convectionForcing[i]);
    }
    const std::array<double,7> sensors{Cast(du),Cast(delta),
        Cast((sample.kN-static_cast<Real>(m_original.K()))/m_original.K()),
        Cast((sample.epsilonN-static_cast<Real>(m_reference.numerator)/sample.heightN)/
             (static_cast<Real>(m_reference.numerator)/sample.heightN)),Cast(fp),Cast(fa),
        Cast(sample.normalVelocityP/us)};
    Real departure=0.0L;
    for (int i=0; i<7; ++i) { departure=std::hypot(departure,std::sqrt(static_cast<Real>(p.weights[i]))*sensors[i]); }
    const double d=Cast(departure);
    auto result=EvaluateFixedWeight(sample.matchingVelocity,forcing,AWFEquilibriumWeight(d,p.on,p.off),policy,highRe);
    result.departure=d; result.sensors=sensors;
    return result;
}
} // namespace Nektar::SolverUtils
