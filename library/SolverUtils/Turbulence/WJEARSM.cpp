///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <SolverUtils/Turbulence/WJEARSM.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace Nektar::SolverUtils
{
namespace
{
using Real = long double;
// First-order directional arithmetic shares every tensor and coefficient
// operation between the value and action. The root itself is solved in reals.
struct Jet
{
    Real v = 0, d = 0;
    Jet() = default;
    Jet(Real value, Real direction = 0) : v(value), d(direction) {}
};
Jet operator+(Jet a, Jet b) { return {a.v+b.v, a.d+b.d}; }
Jet operator-(Jet a, Jet b) { return {a.v-b.v, a.d-b.d}; }
Jet operator-(Jet a) { return {-a.v, -a.d}; }
Jet operator*(Jet a, Jet b) { return {a.v*b.v, a.d*b.v+a.v*b.d}; }
Jet operator/(Jet a, Jet b) { return {a.v/b.v, (a.d-a.v/b.v*b.d)/b.v}; }
using Matrix = std::array<std::array<Jet, 3>, 3>;
Matrix Product(const Matrix &a, const Matrix &b)
{
    Matrix c{};
    for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j)
            for (int l=0; l<3; ++l) c[i][j] = c[i][j]+a[i][l]*b[l][j];
    return c;
}
Jet Trace(const Matrix &a) { return a[0][0]+a[1][1]+a[2][2]; }
void Finite(double v)
{
    if (!std::isfinite(v)) throw std::invalid_argument("Nonfinite WJ input");
}
void Positive(double v)
{
    Finite(v);
    if (!(v>0)) throw std::invalid_argument("WJ requires positive k, epsilon and coefficients");
}
double Cast(Real v)
{
    const double x = static_cast<double>(v);
    if (!std::isfinite(x)) throw std::range_error("Unrepresentable WJ result");
    return x;
}
constexpr Real cp = 9.0L/5.0L;
Real Root(Real s, Real w)
{
    if (!(std::isfinite(s) && std::isfinite(w) && s>=0 && w<=0))
        throw std::range_error("Invalid WJ invariants");
    if (s==0) return cp;
    // F/N is strictly increasing on N>=cp. The w<=0 term bounds the
    // root by the positive root of N(N-cp)=2.7s, irrespective of rotation.
    Real low=cp, high=cp+std::sqrt(2.7L*s);
    auto f = [=](Real n) { return (n-cp)*(n-2*w/n)-2.7L*s; };
    if (!std::isfinite(high)) throw std::range_error("Unbounded WJ root");
    for (int i=0; i<256; ++i)
    {
        const Real mid=low+(high-low)/2;
        if (mid==low || mid==high) break;
        if (f(mid)>0) high=mid; else low=mid;
    }
    return low+(high-low)/2;
}

RANSClosureValues Compute(double kval, double eval, const RANSTensor &g,
    double dk, double de, const RANSTensor &dg, const WJParameters &p,
    WJGradientPolicy policy, double nu, bool direction)
{
    p.Validate(); Positive(kval); Positive(eval); Finite(dk); Finite(de); Finite(nu);
    if (nu!=0) throw std::invalid_argument("WJ high-Re outer diffusion requires outerNu=0");
    if (policy!=WJGradientPolicy::Incompressible && policy!=WJGradientPolicy::DevStrainRawWorkV1)
        throw std::invalid_argument("Unknown WJ gradient completion");
    Jet k{kval,dk}, e{eval,de}, tau=k/e, theta;
    Matrix raw{}, s{}, o{}, dev{};
    Real norm=0, dnorm=0;
    for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j)
        {
            Finite(g[i][j]); Finite(dg[i][j]);
            raw[i][j]={g[i][j],dg[i][j]};
            norm=std::max(norm,std::abs(Real(g[i][j])));
            dnorm=std::max(dnorm,std::abs(Real(dg[i][j])));
        }
    theta=Trace(raw);
    const Real tol=64*std::numeric_limits<double>::epsilon();
    if (policy==WJGradientPolicy::Incompressible &&
        (std::abs(theta.v)>tol*norm || std::abs(theta.d)>tol*dnorm))
        throw std::invalid_argument("Strict WJ core requires trace-free gradient and direction");
    for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j)
        {
            dev[i][j]=raw[i][j]-(i==j ? theta/3 : Jet{});
            s[i][j]=tau*((raw[i][j]+raw[j][i])/2-(i==j ? theta/3 : Jet{}));
            o[i][j]=tau*(raw[i][j]-raw[j][i])/2;
        }
    const auto s2=Product(s,s), o2=Product(o,o);
    const Jet is=Trace(s2), iw=Trace(o2), iv=Trace(Product(s,o2));
    const Jet i5=Trace(Product(s2,o2));
    const Real nr=Root(is.v,iw.v);
    const Real fn=cp*(nr*nr-2*iw.v)/nr+2*nr*(nr-cp);
    const Jet n{nr,(2.7L*nr*is.d+2*(nr-cp)*iw.d)/fn};
    const Jet fnj=cp*(n*n-2*iw)/n+2*n*(n-cp);
    const Jet q=(5.0L/6)*(n*n-2*iw)*(2*n*n-iw);
    const std::array<Jet,5> beta{{-n*(2*n*n-7*iw)/q, -12*iv/(n*q),
                                -2*(n*n-2*iw)/q, -6*n/q, 6/q}};
    const Jet cmu=3*n/(5*(n*n-2*iw)), nut=cmu*k*tau, nud=p.cd*k*tau;
    if (!(Cast(nut.v)>0 && Cast(nud.v)>0)) throw std::range_error("WJ viscosity underflow");
    const auto so=Product(s,o), os=Product(o,s), so2=Product(s,o2), o2s=Product(o2,s);
    const auto oso2=Product(os,o2), o2so=Product(o2s,o);
    Matrix a{}, r{}, remainder{}, stress{};
    Jet production, pdev;
    for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j)
        {
            const Jet delta=i==j ? 1 : 0;
            a[i][j]=beta[0]*s[i][j]+beta[1]*(o2[i][j]-iw*delta/3)+
                beta[2]*(so[i][j]-os[i][j])+
                beta[3]*(so2[i][j]+o2s[i][j]-2*iv*delta/3)+
                beta[4]*(oso2[i][j]-o2so[i][j]);
            stress[i][j]=k*a[i][j];
            r[i][j]=k*(2*delta/3+a[i][j]);
            remainder[i][j]=stress[i][j]+2*nut*s[i][j]/tau;
            production=production-r[i][j]*raw[i][j];
            pdev=pdev-stress[i][j]*dev[i][j];
        }
    const Jet ep=p.cEpsilon1*e/k*production, ed=p.cEpsilon2*e*e/k;
    auto get=[=](Jet x) { return Cast(direction ? x.d : x.v); };
    RANSClosureValues v{};
    v.cMu=get(cmu); v.eddyViscosity=get(nut); v.scalarViscosity=get(nud);
    v.strainInvariant=get(is); v.rotationInvariant=get(iw); v.divergence=get(theta);
    v.production=get(production); v.productionDev=get(pdev); v.kSource=get(production-e);
    v.epsilonProduction=get(ep); v.epsilonDestruction=get(ed); v.epsilonSource=get(ep-ed);
    v.momentumDiffusion=get(nut); v.kDiffusion=get(nud/p.sigmaK); v.epsilonDiffusion=get(nud/p.sigmaEpsilon);
    for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j)
        {
            v.anisotropy[i][j]=get(a[i][j]); v.reynoldsStress[i][j]=get(r[i][j]);
            v.deviatoricStress[i][j]=get(stress[i][j]); v.stressRemainder[i][j]=get(remainder[i][j]);
        }
    v.wj.n=get(n); v.wj.s=get(is); v.wj.w=get(iw); v.wj.i4=get(iv); v.wj.i5=get(i5);
    v.wj.q=get(q); v.wj.rootDerivative=get(fnj);
    v.wj.rootResidual=get((n-cp)*(n*n-2*iw)-2.7L*is*n);
    for (int i=0; i<5; ++i) v.wj.beta[i]=get(beta[i]);
    return v;
}
}

void WJParameters::Validate() const
{
    for (double x : {c1,c2,cd,cEpsilon1,cEpsilon2,sigmaK,sigmaEpsilon}) Positive(x);
    if (c1!=1.8 || c2!=5.0/9.0)
        throw std::invalid_argument("WJ-1 supports only c1=1.8 and c2=5/9");
}
RANSClosureValues EvaluateWJClosure(double k, double e, const RANSTensor &g,
    const WJParameters &p, WJGradientPolicy policy, double nu)
{ return Compute(k,e,g,0,0,{},p,policy,nu,false); }
RANSClosureValues ApplyWJClosureDirection(double k, double e, const RANSTensor &g,
    double dk, double de, const RANSTensor &dg, const WJParameters &p,
    WJGradientPolicy policy, double nu)
{ return Compute(k,e,g,dk,de,dg,p,policy,nu,true); }

RANSVector WJOuterTraction(const RANSTensor &r, double p, const RANSVector &n)
{
    Finite(p);
    Real norm=0;
    for (double x:n) { Finite(x); norm+=Real(x)*x; }
    if (std::abs(norm-1)>64*std::numeric_limits<double>::epsilon())
        throw std::invalid_argument("WJ traction requires a unit normal");
    RANSVector t{};
    for (int i=0; i<3; ++i)
    {
        Real component=-Real(p)*n[i];
        for (int j=0; j<3; ++j) { Finite(r[i][j]); component-=Real(r[i][j])*n[j]; }
        t[i]=Cast(component);
    }
    return t;
}
double WJInterfacePressureJump(double k, double d, const RANSTensor &m,
    const RANSTensor &r, const RANSVector &n)
{
    Finite(k); Positive(d);
    if (k<0) throw std::invalid_argument("WJ wall pressure requires nonnegative local k");
    const auto outer=WJOuterTraction(r,0,n); // Also validates R and n.
    Real jump=-2*Real(k)/3;
    for (int i=0; i<3; ++i)
    {
        jump-=Real(n[i])*outer[i];
        for (int j=0; j<3; ++j) { Finite(m[i][j]); jump+=Real(d)*n[i]*m[i][j]*n[j]; }
    }
    return Cast(jump);
}
double WJMechanicalPressure(double pi, double k, double potential)
{
    Finite(pi); Finite(k); Finite(potential);
    if (k<0) throw std::invalid_argument("WJ pressure recovery requires nonnegative local k");
    return Cast(Real(pi)-2*Real(k)/3-potential);
}
RANSClosureKernel::RANSClosureKernel(const WJParameters &p, WJGradientPolicy policy)
    : m_parameters(p), m_policy(policy) {}
RANSClosureKernel RANSClosureKernel::WJ(const WJParameters &p, WJGradientPolicy policy)
{
    p.Validate();
    if (policy!=WJGradientPolicy::Incompressible && policy!=WJGradientPolicy::DevStrainRawWorkV1)
        throw std::invalid_argument("Unknown WJ factory gradient completion");
    return {p,policy};
}
RANSClosureValues RANSClosureKernel::Evaluate(double k,double e,const RANSTensor &g,double nu) const
{ return EvaluateWJClosure(k,e,g,m_parameters,m_policy,nu); }
RANSClosureValues RANSClosureKernel::Direction(double k,double e,const RANSTensor &g,
    double dk,double de,const RANSTensor &dg,double nu) const
{ return ApplyWJClosureDirection(k,e,g,dk,de,dg,m_parameters,m_policy,nu); }
}
