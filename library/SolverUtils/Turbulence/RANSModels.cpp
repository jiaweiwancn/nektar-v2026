///////////////////////////////////////////////////////////////////////////////
// File: RANSModels.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <SolverUtils/Turbulence/RANSModels.h>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>

namespace Nektar::SolverUtils
{
namespace
{
using Real = long double;
using Matrix = std::array<std::array<Real, 3>, 3>;

void Finite(double value, const char *name)
{
    if (!std::isfinite(value)) { throw std::invalid_argument(name); }
}
void Positive(double value, const char *name)
{
    Finite(value, name);
    if (!(value > 0.0)) { throw std::invalid_argument(name); }
}
double Cast(Real value)
{
    const auto result = static_cast<double>(value);
    if (!std::isfinite(result)) { throw std::range_error("RANS result is not representable"); }
    return result;
}
Real Coefficient(Real q)
{
    const Real factor = q > 8.0L ? 1.0L : -std::expm1(-0.36L * std::exp(0.75L*q));
    if (q > 1.0L)
    {
        const Real reciprocal = (1.0L/q) / std::sqrt(q);
        return 0.3L * factor * reciprocal / (reciprocal + 0.35L);
    }
    return 0.3L * factor / (1.0L + 0.35L*q*std::sqrt(q));
}
Real Yap(Real r, const RANSParameters &p)
{
    return p.yapCoefficient * std::max(0.0L, r-1.0L) * r*r;
}
Matrix Product(const Matrix &a, const Matrix &b)
{
    Matrix c{};
    for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j)
            for (int k=0; k<3; ++k) { c[i][j] += a[i][k]*b[k][j]; }
    return c;
}
Real Root(const std::function<Real(Real)> &function)
{
    Real low=0.0L, high=1.0L;
    while (function(high) < 0.0L)
    {
        high *= 2.0L;
        if (!std::isfinite(high)) { throw std::range_error("RANS reference root is unbounded"); }
    }
    for (int iteration=0; iteration<256; ++iteration)
    {
        const Real middle=low+(high-low)/2.0L;
        if (middle==low || middle==high) { break; }
        const Real value=function(middle);
        if (!std::isfinite(value)) { throw std::range_error("Invalid RANS root function"); }
        if (value > 0.0L) { high=middle; } else { low=middle; }
    }
    return low+(high-low)/2.0L;
}
} // namespace

void RANSParameters::Validate() const
{
    for (double value : {cMu,cEpsilon1,cEpsilon2,sigmaK,sigmaEpsilon,yapLength})
        { Positive(value, "RANS coefficients must be finite and positive"); }
    Finite(yapCoefficient, "Invalid Yap coefficient");
    if (yapCoefficient < 0.0) { throw std::invalid_argument("Negative Yap coefficient"); }
}

double CLSCmu(double q)
{
    Finite(q, "Nonfinite CLS invariant");
    if (q < 0.0) { throw std::invalid_argument("Negative CLS invariant"); }
    return Cast(Coefficient(q));
}

double CLSYapFunction(double ratio, const RANSParameters &p)
{
    p.Validate();
    Finite(ratio, "Nonfinite Yap length ratio");
    if (ratio < 0.0) { throw std::invalid_argument("Negative Yap length ratio"); }
    return Cast(Yap(ratio,p));
}

std::array<double, 7> CLSCubicCoefficients(double cmu)
{
    Positive(cmu, "Invalid cubic Cmu");
    const Real square=static_cast<Real>(cmu)*cmu;
    return {-0.1,0.1,0.26,Cast(-10.0L*square),0.0,
            Cast(-5.0L*square),Cast(5.0L*square)};
}

RANSClosureValues EvaluateRANSClosure(
    RANSModel model, double k, double epsilon, const RANSTensor &gradient,
    double distance, const RANSParameters &p, double nu)
{
    p.Validate();
    Positive(k, "RANS k must be finite and positive");
    Positive(epsilon, "RANS dissipation must be finite and positive");
    Finite(nu, "Nonfinite resolved viscosity");
    if (nu < 0.0) { throw std::invalid_argument("Negative resolved viscosity"); }
    if (model!=RANSModel::StandardKEpsilon && model!=RANSModel::CLSHighRe)
        { throw std::invalid_argument("Unknown RANS closure"); }
    if (model==RANSModel::CLSHighRe) { Positive(distance, "CLS requires geometric wall distance > 0"); }

    const Real time=static_cast<Real>(k)/epsilon;
    Matrix s{},o{};
    Real ss=0.0L,oo=0.0L,divergence=0.0L;
    for (int i=0; i<3; ++i)
    {
        divergence+=gradient[i][i];
        for (int j=0; j<3; ++j)
        {
            Finite(gradient[i][j], "Nonfinite velocity gradient");
            s[i][j]=time*(static_cast<Real>(gradient[i][j])+gradient[j][i]);
            o[i][j]=time*(static_cast<Real>(gradient[i][j])-gradient[j][i]);
            ss+=s[i][j]*s[i][j]; oo+=o[i][j]*o[i][j];
        }
    }
    const Real strain=std::sqrt(ss/2.0L),rotation=std::sqrt(oo/2.0L);
    const Real cmu=model==RANSModel::CLSHighRe ? Coefficient(std::max(strain,rotation)) : p.cMu;
    const Real nut=cmu*k*time;
    RANSClosureValues result{};
    result.cMu=Cast(cmu); result.eddyViscosity=Cast(nut);
    if (!(result.eddyViscosity>0.0)) { throw std::range_error("Eddy viscosity underflow"); }
    result.strainInvariant=Cast(strain); result.rotationInvariant=Cast(rotation);
    result.divergence=Cast(divergence);
    Matrix quadratic{},cubic{};
    if (model==RANSModel::CLSHighRe)
    {
        const Matrix s2=Product(s,s),o2=Product(o,o),os=Product(o,s),so=Product(s,o);
        const Matrix s2o=Product(s2,o),os2=Product(o,s2);
        const Real c4=-10.0L*cmu*cmu,c6=-5.0L*cmu*cmu,c7=5.0L*cmu*cmu;
        for (int i=0; i<3; ++i)
            for (int j=0; j<3; ++j)
            {
                const Real delta=i==j ? 1.0L : 0.0L;
                // O O^T = -O^2; published c5 is exactly zero.
                quadratic[i][j]=k*cmu*(-0.1L*(s2[i][j]-ss*delta/3.0L)+
                    0.1L*(os[i][j]-so[i][j])+0.26L*(-o2[i][j]-oo*delta/3.0L));
                cubic[i][j]=k*cmu*(c4*(s2o[i][j]-os2[i][j])+
                                               (c6*ss+c7*oo)*s[i][j]);
            }
    }
    Real production=0.0L;
    for (int i=0; i<3; ++i)
        for (int j=0; j<3; ++j)
        {
            const Real stress=(i==j ? 2.0L*k/3.0L : 0.0L)-k*cmu*s[i][j]+
                              quadratic[i][j]+cubic[i][j];
            result.reynoldsStress[i][j]=Cast(stress);
            result.quadraticStress[i][j]=Cast(quadratic[i][j]);
            result.cubicStress[i][j]=Cast(cubic[i][j]);
            production-=stress*gradient[i][j];
        }
    const Real epsScale=static_cast<Real>(epsilon)*epsilon/k;
    const Real ratio=model==RANSModel::CLSHighRe ?
        static_cast<Real>(k)*std::sqrt(static_cast<Real>(k))/(p.yapLength*epsilon*static_cast<Real>(distance)) : 0.0L;
    const Real y=model==RANSModel::CLSHighRe ? epsScale*Yap(ratio,p) : 0.0L;
    result.production=Cast(production);
    result.kSource=Cast(production-epsilon);
    result.epsilonProduction=Cast(static_cast<Real>(p.cEpsilon1)*epsilon/k*production);
    result.epsilonDestruction=Cast(p.cEpsilon2*epsScale);
    result.yap=Cast(y);
    result.epsilonSource=Cast(static_cast<Real>(p.cEpsilon1)*epsilon/k*production-p.cEpsilon2*epsScale+y);
    result.momentumDiffusion=Cast(nu+nut);
    result.kDiffusion=Cast(nu+nut/p.sigmaK);
    result.epsilonDiffusion=Cast(nu+nut/p.sigmaEpsilon);
    return result;
}

double RANSTotalDissipation(double k,double modified,double nu,const RANSVector &grad)
{
    Positive(k,"Invalid k for dissipation split"); Positive(modified,"Invalid modified dissipation");
    Finite(nu,"Nonfinite molecular viscosity");
    if (nu<0.0) { throw std::invalid_argument("Negative molecular viscosity"); }
    Real norm=0.0L;
    for (double value : grad) { Finite(value,"Nonfinite k gradient"); norm+=static_cast<Real>(value)*value; }
    return Cast(modified+static_cast<Real>(nu)*norm/(2.0L*k));
}

RANSNeutralReference::RANSNeutralReference(
    RANSModel model,const RANSParameters &p,double q,double kap,double cmu)
    : m_model(model),m_parameters(p),m_q(q),m_kappa(kap),m_cmu(cmu) {}

RANSNeutralReference RANSNeutralReference::Standard(const RANSParameters &p)
{
    p.Validate();
    if (!(p.cEpsilon2>p.cEpsilon1)) { throw std::invalid_argument("Neutral reference requires Ceps2 > Ceps1"); }
    const Real q=1.0L/std::sqrt(static_cast<Real>(p.cMu));
    return {RANSModel::StandardKEpsilon,p,Cast(q),
            Cast(std::sqrt((p.cEpsilon2-p.cEpsilon1)*p.sigmaEpsilon/q)),p.cMu};
}

RANSNeutralReference RANSNeutralReference::CLS(const RANSParameters &p)
{
    p.Validate();
    if (!(p.cEpsilon2>p.cEpsilon1)) { throw std::invalid_argument("Neutral reference requires Ceps2 > Ceps1"); }
    const Real q=Root([](Real x){ return Coefficient(x)*x*x-1.0L; });
    const Real kap=Root([&](Real x) {
        return q*x*x/p.sigmaEpsilon+Yap(x*q*std::sqrt(q)/p.yapLength,p)-
               (p.cEpsilon2-p.cEpsilon1);
    });
    return {RANSModel::CLSHighRe,p,Cast(q),Cast(kap),Cast(Coefficient(q))};
}

NeutralABLValues RANSNeutralReference::AtHeight(double z,double u,double H,double U) const
{
    Positive(z,"Invalid reference height"); Positive(u,"Invalid friction velocity");
    Positive(H,"Invalid matching height"); Finite(U,"Invalid matching velocity");
    const Real star=u,eps=star*star*star/(m_kappa*static_cast<Real>(z));
    return {Cast(U+star/m_kappa*(std::log(static_cast<Real>(z))-std::log(static_cast<Real>(H)))),
            Cast(m_q*star*star),Cast(eps),Cast(m_kappa*star*z),
            Cast(star/(m_kappa*static_cast<Real>(z))),Cast(-eps/z),
            Cast(-star*star*star*star/(m_parameters.sigmaEpsilon*static_cast<Real>(z))),
            Cast(star*star)};
}

NeutralABLValues RANSNeutralReference::StandardShifted(double height,double u,double z0) const
{
    if (m_model!=RANSModel::StandardKEpsilon) { throw std::invalid_argument("CLS uses geometric wall distance without a roughness shift"); }
    Finite(height,"Invalid physical height"); Positive(z0,"Invalid aerodynamic roughness");
    if (height<0.0) { throw std::invalid_argument("Negative physical height"); }
    return AtHeight(Cast(static_cast<Real>(height)+z0),u,z0,0.0);
}

std::array<double,2> RANSNeutralReference::StandardLogWallTraction(
    double k,const std::array<double,2> &u,double height,double z0,double shift) const
{
    if (m_model!=RANSModel::StandardKEpsilon) { throw std::invalid_argument("Standard log wall requires standard k-epsilon"); }
    Positive(k,"Invalid wall k"); Positive(z0,"Invalid aerodynamic roughness");
    Finite(height,"Invalid wall sample height"); Finite(shift,"Invalid height shift");
    const Real z=static_cast<Real>(height)+shift;
    if (!(z>z0)) { throw std::invalid_argument("Log-wall matching height must exceed its intercept"); }
    const Real factor=std::pow(static_cast<Real>(m_cmu),0.25L)*std::sqrt(static_cast<Real>(k))*m_kappa/
                      (std::log(z)-std::log(static_cast<Real>(z0)));
    for (double value : u) { Finite(value,"Invalid wall velocity"); }
    return {Cast(factor*u[0]),Cast(factor*u[1])};
}
} // namespace Nektar::SolverUtils
