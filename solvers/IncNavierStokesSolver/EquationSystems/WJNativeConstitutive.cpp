// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeConstitutive.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>

namespace Nektar::WJNative
{
namespace
{
using namespace SolverUtils;
constexpr double Cp=1.8;
struct Dual
{
    double v=0,d=0;
    Dual()=default;
    Dual(double value,double direction=0):v(value),d(direction){}
};
inline Dual operator+(Dual a,Dual b){return {a.v+b.v,a.d+b.d};}
inline Dual operator-(Dual a,Dual b){return {a.v-b.v,a.d-b.d};}
inline Dual operator-(Dual a){return {-a.v,-a.d};}
inline Dual operator*(Dual a,Dual b){return {a.v*b.v,a.d*b.v+a.v*b.d};}
inline Dual operator/(Dual a,Dual b){const double q=a.v/b.v;return {q,(a.d-q*b.d)/b.v};}
inline double Value(double x){return x;}
inline double Value(Dual x){return x.v;}
inline double Direction(double){return 0;}
inline double Direction(Dual x){return x.d;}
template<class T>T Make(double value,double direction)
{if constexpr(std::is_same_v<T,Dual>)return {value,direction};else return value;}
template<class T>using Matrix=std::array<std::array<T,3>,3>;
template<class T>Matrix<T> Product(const Matrix<T>&a,const Matrix<T>&b)
{
    Matrix<T> c{};
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)
        c[i][j]=a[i][0]*b[0][j]+a[i][1]*b[1][j]+a[i][2]*b[2][j];
    return c;
}
template<class T>T Trace(const Matrix<T>&a){return a[0][0]+a[1][1]+a[2][2];}
void Finite(double value)
{if(!std::isfinite(value))throw std::invalid_argument("Nonfinite fast WJ input");}
void Validate(double k,double e,const RANSTensor&g,double dk,double de,const RANSTensor&dg,const WJParameters&p)
{
    p.Validate();Finite(k);Finite(e);Finite(dk);Finite(de);
    if(!(k>0&&e>0))throw std::invalid_argument("Fast WJ requires positive k and epsilon");
    for(int i=0;i<3;++i)for(int j=0;j<3;++j){Finite(g[i][j]);Finite(dg[i][j]);}
}
bool Moderate(double x,double low=1e-80,double high=1e80)
{return x==0||(std::isfinite(x)&&std::abs(x)>=low&&std::abs(x)<=high);}
bool FastScale(double k,double e,const RANSTensor&g,double dk,double de,const RANSTensor&dg,const WJParameters&p)
{
    if(!Moderate(k,1e-100,1e100)||!Moderate(e,1e-100,1e100)||
       !Moderate(k/e,1e-50,1e50)||!Moderate(dk)||!Moderate(de))return false;
    for(double c:{p.cd,p.cEpsilon1,p.cEpsilon2,p.sigmaK,p.sigmaEpsilon})
        if(!Moderate(c,1e-40,1e40))return false;
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)if(!Moderate(g[i][j])||!Moderate(dg[i][j]))return false;
    return true;
}
void RootDiagnostics(double n,double s,double w,double derivative,FastWJStatus&status)
{
    using L=long double;
    const L left=(L(n)-9.L/5)*(L(n)*n-2*L(w)),right=2.7L*s*n;
    const L scale=std::abs(L(n)*n*n)+(9.L/5)*n*n+
        2*std::abs(L(w))*(std::abs(L(n))+9.L/5)+std::abs(right);
    status.normalizedRootResidual=double(std::abs(left-right)/scale);
    status.relativeRootCondition=double(L(derivative)/(L(n)*n+2*std::abs(L(w))+2.7L*s));
}
bool Root(double s,double w,double&n,FastWJStatus&status)
{
    if(!(std::isfinite(s)&&std::isfinite(w)&&s>=0&&w<=0))
    {status.fallback=FastWJFallback::InvariantRange;return false;}
    if(s==0){n=Cp;RootDiagnostics(n,s,w,Cp*(n*n-2*w)/n,status);return true;}
    double low=Cp,high=Cp+std::sqrt(2.7*s);n=high;
    if(!std::isfinite(high)){status.fallback=FastWJFallback::UnverifiedRoot;return false;}
    // F/N is strictly increasing for N>=Cp. Safeguarded Newton retains its
    // physical bracket; the finite fast-path work cap delegates to old Root.
    for(int iteration=0;iteration<20;++iteration)
    {
        ++status.rootIterations;
        const double f=(n-Cp)*(n-2*w/n)-2.7*s;
        const double derivative=2*n-Cp-2*Cp*w/(n*n);
        if(!(std::isfinite(f)&&std::isfinite(derivative)&&derivative>0))break;
        if(std::abs(f)/derivative<=4*std::numeric_limits<double>::epsilon()*n)
        {
            const double fn=Cp*(n*n-2*w)/n+2*n*(n-Cp);
            RootDiagnostics(n,s,w,fn,status);
            if(n>=Cp&&status.normalizedRootResidual<=64*std::numeric_limits<double>::epsilon()&&
               std::isfinite(status.relativeRootCondition)&&status.relativeRootCondition>128*std::numeric_limits<double>::epsilon())return true;
            break;
        }
        if(f>0)high=n;else low=n;
        double next=n-f/derivative;
        if(!std::isfinite(next)||next<=low||next>=high)next=low+(high-low)/2;
        if(next==n)break;
        n=next;
    }
    status.fallback=FastWJFallback::UnverifiedRoot;return false;
}
bool FiniteOutputs(const RANSClosureValues&v)
{
    for(double x:{v.cMu,v.eddyViscosity,v.scalarViscosity,v.strainInvariant,v.rotationInvariant,
        v.divergence,v.production,v.productionDev,v.kSource,v.epsilonProduction,v.epsilonDestruction,
        v.epsilonSource,v.momentumDiffusion,v.kDiffusion,v.epsilonDiffusion,v.wj.n,v.wj.s,v.wj.w,
        v.wj.i4,v.wj.i5,v.wj.q,v.wj.rootDerivative,v.wj.rootResidual})if(!std::isfinite(x))return false;
    for(double x:v.wj.beta)if(!std::isfinite(x))return false;
    for(const auto*t:{&v.anisotropy,&v.reynoldsStress,&v.deviatoricStress,&v.stressRemainder})
        for(const auto&row:*t)for(double x:row)if(!std::isfinite(x))return false;
    return true;
}
template<class T>bool Compute(double kval,double eval,const RANSTensor&g,double dk,double de,
    const RANSTensor&dg,const WJParameters&p,FastWJPair&out)
{
    const T k=Make<T>(kval,dk),e=Make<T>(eval,de),tau=k/e;
    Matrix<T>raw{},s{},o{},dev{};
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)raw[i][j]=Make<T>(g[i][j],dg[i][j]);
    const T theta=Trace(raw);
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)
    {
        dev[i][j]=raw[i][j]-(i==j?theta/T(3):T(0));
        s[i][j]=tau*((raw[i][j]+raw[j][i])/T(2)-(i==j?theta/T(3):T(0)));
        o[i][j]=tau*(raw[i][j]-raw[j][i])/T(2);
        for(T x:{s[i][j],o[i][j]})
            if(!Moderate(Value(x),1e-40,1e20)||!Moderate(Direction(x),1e-40,1e20))
            {out.status.fallback=FastWJFallback::ExtremeScale;return false;}
    }
    const auto s2=Product(s,s),o2=Product(o,o);
    const T is=Trace(s2),iw=Trace(o2),iv=Trace(Product(s,o2)),i5=Trace(Product(s2,o2));
    double nr=0;if(!Root(Value(is),Value(iw),nr,out.status))return false;
    const double fn=Cp*(nr*nr-2*Value(iw))/nr+2*nr*(nr-Cp);
    const T n=Make<T>(nr,(2.7*nr*Direction(is)+2*(nr-Cp)*Direction(iw))/fn);
    const T fnj=T(Cp)*(n*n-T(2)*iw)/n+T(2)*n*(n-T(Cp));
    const T q=T(5./6)*(n*n-T(2)*iw)*(T(2)*n*n-iw);
    const std::array<T,5> beta{{-n*(T(2)*n*n-T(7)*iw)/q,-T(12)*iv/(n*q),
                               -T(2)*(n*n-T(2)*iw)/q,-T(6)*n/q,T(6)/q}};
    // Preserve the exact selected stress Cmu; Cd is independently supplied.
    const T cmu=T(3)*n/(T(5)*(n*n-T(2)*iw)),nut=cmu*k*tau,nud=T(p.cd)*k*tau;
    const auto so=Product(s,o),os=Product(o,s),so2=Product(s,o2),o2s=Product(o2,s);
    const auto oso2=Product(os,o2),o2so=Product(o2s,o);
    Matrix<T> anisotropy{},stress{},covariance{},remainder{};T production{},pdev{};
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)
    {
        const T delta=i==j?T(1):T(0);
        anisotropy[i][j]=beta[0]*s[i][j]+beta[1]*(o2[i][j]-iw*delta/T(3))+
            beta[2]*(so[i][j]-os[i][j])+beta[3]*(so2[i][j]+o2s[i][j]-T(2)*iv*delta/T(3))+
            beta[4]*(oso2[i][j]-o2so[i][j]);
        stress[i][j]=k*anisotropy[i][j];covariance[i][j]=k*(T(2)*delta/T(3)+anisotropy[i][j]);
        remainder[i][j]=stress[i][j]+T(2)*nut*s[i][j]/tau;
        production=production-covariance[i][j]*raw[i][j];pdev=pdev-stress[i][j]*dev[i][j];
    }
    const T ep=T(p.cEpsilon1)*e/k*production,ed=T(p.cEpsilon2)*e*e/k;
    auto output=[&](bool tangent)
    {
        auto get=[&](T x){return tangent?Direction(x):Value(x);};RANSClosureValues v{};
        v.cMu=get(cmu);v.eddyViscosity=get(nut);v.scalarViscosity=get(nud);
        v.strainInvariant=get(is);v.rotationInvariant=get(iw);v.divergence=get(theta);
        v.production=get(production);v.productionDev=get(pdev);v.kSource=get(production-e);
        v.epsilonProduction=get(ep);v.epsilonDestruction=get(ed);v.epsilonSource=get(ep-ed);
        v.momentumDiffusion=get(nut);v.kDiffusion=get(nud/T(p.sigmaK));v.epsilonDiffusion=get(nud/T(p.sigmaEpsilon));
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)
        {v.anisotropy[i][j]=get(anisotropy[i][j]);v.reynoldsStress[i][j]=get(covariance[i][j]);
         v.deviatoricStress[i][j]=get(stress[i][j]);v.stressRemainder[i][j]=get(remainder[i][j]);}
        v.wj.n=get(n);v.wj.s=get(is);v.wj.w=get(iw);v.wj.i4=get(iv);v.wj.i5=get(i5);v.wj.q=get(q);
        v.wj.rootDerivative=get(fnj);v.wj.rootResidual=get((n-T(Cp))*(n*n-T(2)*iw)-T(2.7)*is*n);
        for(int i=0;i<5;++i)v.wj.beta[i]=get(beta[i]);return v;
    };
    out.value=output(false);
    if constexpr(std::is_same_v<T,Dual>)out.direction=output(true);
    if(!FiniteOutputs(out.value)||!FiniteOutputs(out.direction))
    {out.status.fallback=FastWJFallback::NonfiniteOutput;return false;}
    if(!(out.value.eddyViscosity>0&&out.value.scalarViscosity>0))
    {out.status.fallback=FastWJFallback::ViscosityUnderflow;return false;}
    return true;
}
FastWJPair Fallback(double k,double e,const RANSTensor&g,double dk,double de,const RANSTensor&dg,
    const WJParameters&p,bool direction,FastWJStatus status)
{
    FastWJPair result;result.status=status;
    result.value=EvaluateWJClosure(k,e,g,p,WJGradientPolicy::DevStrainRawWorkV1,0);
    if(direction)result.direction=ApplyWJClosureDirection(k,e,g,dk,de,dg,p,WJGradientPolicy::DevStrainRawWorkV1,0);
    RootDiagnostics(result.value.wj.n,result.value.wj.s,result.value.wj.w,result.value.wj.rootDerivative,result.status);
    return result;
}
template<class T>FastWJPair Dispatch(double k,double e,const RANSTensor&g,double dk,double de,
    const RANSTensor&dg,const WJParameters&p)
{
    Validate(k,e,g,dk,de,dg,p);FastWJPair result;
    if(!FastScale(k,e,g,dk,de,dg,p))result.status.fallback=FastWJFallback::ExtremeScale;
    else if(Compute<T>(k,e,g,dk,de,dg,p,result))return result;
    return Fallback(k,e,g,dk,de,dg,p,std::is_same_v<T,Dual>,result.status);
}
}

SolverUtils::RANSClosureValues EvaluateFastWJClosure(double k,double e,const SolverUtils::RANSTensor&g,
    const SolverUtils::WJParameters&p,FastWJStatus*status)
{auto result=Dispatch<double>(k,e,g,0,0,{},p);if(status)*status=result.status;return result.value;}
SolverUtils::RANSClosureValues ApplyFastWJClosureDirection(double k,double e,const SolverUtils::RANSTensor&g,
    double dk,double de,const SolverUtils::RANSTensor&dg,const SolverUtils::WJParameters&p,FastWJStatus*status)
{auto result=Dispatch<Dual>(k,e,g,dk,de,dg,p);if(status)*status=result.status;return result.direction;}
FastWJPair EvaluateFastWJValueAndDirection(double k,double e,const SolverUtils::RANSTensor&g,
    double dk,double de,const SolverUtils::RANSTensor&dg,const SolverUtils::WJParameters&p)
{return Dispatch<Dual>(k,e,g,dk,de,dg,p);}
const char*FastWJFallbackName(FastWJFallback reason)
{
    switch(reason)
    {
        case FastWJFallback::None:return "none";
        case FastWJFallback::ExtremeScale:return "extreme_scale";
        case FastWJFallback::InvariantRange:return "invariant_range";
        case FastWJFallback::UnverifiedRoot:return "unverified_root";
        case FastWJFallback::NonfiniteOutput:return "nonfinite_output";
        case FastWJFallback::ViscosityUnderflow:return "viscosity_underflow";
    }
    return "unknown";
}
}
