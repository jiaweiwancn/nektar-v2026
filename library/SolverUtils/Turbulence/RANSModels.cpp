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
#include <sstream>
#include <iomanip>

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

bool IsNewRANSModel(RANSModel m)
{ return m==RANSModel::A93 || m==RANSModel::A95 || m==RANSModel::E94E95; }
const char *RANSModelName(RANSModel m)
{
    switch(m)
    {
        case RANSModel::StandardKEpsilon:return "StandardKEpsilon";
        case RANSModel::CLSHighRe:return "CLSHighRe";
        case RANSModel::A93:return "A93";
        case RANSModel::A95:return "A95";
        case RANSModel::E94E95:return "E94E95";
    }
    throw std::invalid_argument("Unknown RANS model identity");
}
RANSModel ParseRANSModel(const std::string &name)
{
    for(auto m:{RANSModel::StandardKEpsilon,RANSModel::CLSHighRe,RANSModel::A93,RANSModel::A95,RANSModel::E94E95})
        if(name==RANSModelName(m))return m;
    throw std::invalid_argument("Unknown RANS model identity: "+name);
}
namespace
{
// Symmetric Jacobi eigenvalue evaluation. Avoid acos-domain clipping and any
// coefficient continuation at repeated eigenvalues. Forward values only.
Real LargestStrainEigenvalue(Matrix a)
{
    for(int iteration=0;iteration<64;++iteration)
    {
        int p=0,q=1;
        for(int i=0;i<3;++i)for(int j=i+1;j<3;++j)
            if(std::abs(a[i][j])>std::abs(a[p][q])){p=i;q=j;}
        if(a[p][q]==0)break;
        const Real angle=.5L*std::atan2(2*a[p][q],a[q][q]-a[p][p]);
        const Real c=std::cos(angle),v=std::sin(angle);
        const Real pp=c*c*a[p][p]-2*c*v*a[p][q]+v*v*a[q][q];
        const Real qq=v*v*a[p][p]+2*c*v*a[p][q]+c*c*a[q][q];
        for(int i=0;i<3;++i)if(i!=p && i!=q)
        {const Real x=a[i][p],y=a[i][q];a[i][p]=a[p][i]=c*x-v*y;a[i][q]=a[q][i]=v*x+c*y;}
        a[p][p]=pp;a[q][q]=qq;a[p][q]=a[q][p]=0;
    }
    return std::max({a[0][0],a[1][1],a[2][2]});
}
RANSClosureValues NewClosure(RANSModel model,double k,double epsilon,
    const RANSTensor &gradient,const RANSParameters &p,double nu)
{
    Positive(k,"RANS k must be finite and positive");
    Positive(epsilon,"RANS total epsilon must be finite and positive");p.Validate();
    const bool e=model==RANSModel::E94E95;
    if(nu!=0 || p.yapCoefficient!=0 || p.sigmaK!=1 ||
       p.cEpsilon1!=(e?.43:1.44) || p.cEpsilon2!=(e?1.9:1.92))
        throw std::invalid_argument("New closure requires resolved original high-Re coefficients and total epsilon");
    const Real time=Real(k)/epsilon;
    Matrix s{},o{};Real div=0,ss=0,oo=0;
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)Finite(gradient[i][j],"Nonfinite velocity gradient");
    for(int i=0;i<3;++i)div+=gradient[i][i];
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)
    {
        // Unlike inherited CLS, these are tau*S_d and tau*Omega, without
        // doubling. The isotropic raw-divergence production is retained below.
        s[i][j]=time*((Real(gradient[i][j])+gradient[j][i])/2-(i==j?div/3:0));
        o[i][j]=time*(Real(gradient[i][j])-gradient[j][i])/2;
        ss+=s[i][j]*s[i][j];oo+=o[i][j]*o[i][j];
    }
    const Real eta=std::sqrt(2*ss),xi=std::sqrt(2*oo);
    Real cmu,cq=0;
    if(model==RANSModel::A93)cmu=(2.L/3)/(1.25L+eta+.9L*xi);
    else
    {
        if(ss==0 && oo!=0)
            throw std::domain_error("A95/E94E95 undefined scalar-viscosity limit at zero strain and nonzero rotation");
        const Real invariant=ss==0?0:3*LargestStrainEigenvalue(s)*std::sqrt((ss+oo)/ss);
        cmu=1/((e?4.L:6.5L)+invariant);
        if(!e)
        {
            const Real radicand=1-9*cmu*cmu*ss;
            if(radicand<0)throw std::domain_error("A95 negative square-root radicand; no constitutive clipping");
            // radicand==0 is a forward value, not a smooth tangent promise.
            cq=std::sqrt(radicand)/(1+6*std::sqrt(ss)*std::sqrt(oo));
        }
    }
    RANSClosureValues result{};const Real nut=cmu*k*time;
    result.cMu=Cast(cmu);result.eddyViscosity=Cast(nut);
    result.strainInvariant=Cast(eta);result.rotationInvariant=Cast(xi);result.divergence=Cast(div);
    const auto s2=Product(s,s),o2=Product(o,o),os=Product(o,s),so=Product(s,o);
    Real production=0;
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)
    {
        const Real delta=i==j?1:0;
        const Real q=model==RANSModel::A93?k/(1000+eta*eta*eta)*
            (3*(s2[i][j]-ss*delta/3)-19*(o2[i][j]+oo*delta/3)+15*(os[i][j]-so[i][j])):
            model==RANSModel::A95?2*cq*k*(os[i][j]-so[i][j]):0;
        const Real r=2*Real(k)*delta/3-2*k*cmu*s[i][j]+q;
        result.reynoldsStress[i][j]=Cast(r);result.quadraticStress[i][j]=Cast(q);
        production-=r*gradient[i][j];
    }
    const Real scale=Real(epsilon)*epsilon/k;
    const Real ep=e?std::max(.43L,eta/(eta+5))*eta*scale:Real(p.cEpsilon1)*epsilon/k*production;
    const Real ed=Real(p.cEpsilon2)*scale;
    result.production=Cast(production);result.kSource=Cast(production-epsilon);
    result.epsilonProduction=Cast(ep);result.epsilonDestruction=Cast(ed);result.epsilonSource=Cast(ep-ed);
    result.momentumDiffusion=Cast(nut);result.kDiffusion=Cast(nut/p.sigmaK);result.epsilonDiffusion=Cast(nut/p.sigmaEpsilon);
    return result;
}
}

static RANSClosureValues EvaluateRANSClosureSelected(
    RANSModel model, double k, double epsilon, const RANSTensor &gradient,
    double distance, const RANSParameters &p, double nu, int branch)
{
    if(IsNewRANSModel(model))
    {
        if(branch!=0)throw std::invalid_argument("A93/A95/E94E95 constitutive tangent action unsupported; forward one-sweep path only");
        return NewClosure(model,k,epsilon,gradient,p,nu);
    }
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
    const Real cmu=model==RANSModel::CLSHighRe ? Coefficient(branch == 1 ? strain : branch == 2 ? rotation : std::max(strain,rotation)) : p.cMu;
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

RANSClosureValues EvaluateRANSClosure(
    RANSModel model,double k,double epsilon,const RANSTensor &gradient,
    double distance,const RANSParameters &p,double nu)
{
    return EvaluateRANSClosureSelected(model,k,epsilon,gradient,distance,p,nu,0);
}

RANSClosureValues EvaluateRANSClosureLinearization(
    RANSModel model,double k,double epsilon,const RANSTensor &gradient,
    double distance,const RANSParameters &p,double nu,bool strainBranch)
{
    return EvaluateRANSClosureSelected(model,k,epsilon,gradient,distance,p,nu,strainBranch ? 1 : 2);
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

RANSNeutralReference RANSNeutralReference::Resolve(RANSModel model,RANSKappaPolicy policy,
    double prescribed,RANSWallFormulation wall,double molecularNu)
{
    if(!IsNewRANSModel(model))throw std::invalid_argument("Resolve requires A93, A95 or E94E95; retain legacy factories for legacy models");
    if(policy!=RANSKappaPolicy::Natural && policy!=RANSKappaPolicy::Prescribed)
        throw std::invalid_argument("Unknown kappa policy");
    if(wall!=RANSWallFormulation::EPAWF && wall!=RANSWallFormulation::AcceptedBudgetAWF)
        throw std::invalid_argument("Unknown wall formulation");
    Positive(molecularNu,"Physical inner molecular viscosity must be positive");
    const bool e=model==RANSModel::E94E95;
    const Real q=model==RANSModel::A93?Root([](Real q){return (2.L/3)*q*q-1.9L*q-1.25L;}):
        Root([&](Real q){return q*q-3/std::sqrt(2.L)*q-(e?4.L:6.5L);});
    const Real delta=e?1.9L-.43L*q:.48L;
    RANSParameters p;p.yapCoefficient=0;p.sigmaK=1;p.cEpsilon1=e?.43:1.44;p.cEpsilon2=e?1.9:1.92;
    p.sigmaEpsilon=e?1.2:1.3;
    Real kap=std::sqrt(Real(p.sigmaEpsilon)*delta/q);
    if(policy==RANSKappaPolicy::Prescribed)
    {Positive(prescribed,"Prescribed kappa must be positive");kap=prescribed;p.sigmaEpsilon=Cast(kap*kap*q/delta);}
    RANSNeutralReference result(model,p,Cast(q),Cast(kap),Cast(1/(q*q)));
    result.m_policy=policy;result.m_wall=wall;result.m_molecularNu=molecularNu;return result;
}
std::string RANSNeutralReference::Signature() const
{
    std::ostringstream s;s<<std::setprecision(17)<<RANSModelName(m_model)<<":OriginalHighRe-v1";
    s<<":q="<<m_q<<":kappa="<<m_kappa<<":policy="<<(m_policy==RANSKappaPolicy::Natural?"Natural":"Prescribed")
     <<":sigmaK="<<m_parameters.sigmaK<<":sigmaE="<<m_parameters.sigmaEpsilon
     <<":Ceps1="<<m_parameters.cEpsilon1<<":Ceps2="<<m_parameters.cEpsilon2
     <<":wall="<<(m_wall==RANSWallFormulation::EPAWF?"EPAWF:SCI-original-.09-2.55":"AcceptedBudgetAWF:compatible-q-kappa-v1")
     <<":outerNu=0:innerNu="<<m_molecularNu<<":eps_tilde=TotalEpsilon:DeviatoricStrainV1"
     <<":NoRegularization:NoYap:source="<<(m_model==RANSModel::E94E95?"C1(eta)*S*epsilon":"Ceps1*epsilon*P/k")
     <<":wallSink="<<(m_wall==RANSWallFormulation::EPAWF?"P-Dpre:redistribution-once":"P-Ew(K):frozen-prescribed-zero-reaction");
    return s.str();
}
RANSClosureValues EvaluateRANSClosure(const RANSResolvedModel &model,double k,double epsilon,
    const RANSTensor &g,double distance)
{return EvaluateRANSClosure(model.Model(),k,epsilon,g,distance,model.Parameters(),0);}

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
