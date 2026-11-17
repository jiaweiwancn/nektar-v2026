// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeWallV2.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeConstitutive.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace Nektar::WJNative
{
namespace
{
void Require(bool value,const char*text){if(!value)throw std::runtime_error(text);}
struct Data{Tensor strain{},iso{},projector{};Vector mn{};double chi=0,dchi=0;};
Data Inputs(double k,double e,const Tensor&g,const Vector&n,double nu,double a,double r,double H)
{
    Require(std::isfinite(k)&&k>0&&std::isfinite(e)&&e>0&&std::isfinite(nu)&&nu>0&&std::isfinite(a)&&a>=0&&std::isfinite(r)&&r>=0&&std::isfinite(H)&&H>0,"invalid v2 wall inputs");
    double norm=0;for(double x:n){Require(std::isfinite(x),"nonfinite wall normal");norm+=x*x;}Require(std::abs(norm-1)<1e-12,"v2 wall normal is not unit");Data d;
    for(int i=0;i<3;++i)for(int j=0;j<3;++j){Require(std::isfinite(g[i][j]),"nonfinite wall gradient");d.projector[i][j]=(i==j?1.:0.)-n[i]*n[j];d.strain[i][j]=g[i][j]+g[j][i];d.mn[i]+=d.strain[i][j]*n[j];}
    Vector b{};for(int i=0;i<3;++i)for(int j=0;j<3;++j)b[i]-=a*d.projector[i][j]*d.mn[j];
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)d.iso[i][j]=(i==j?2*k/3:0)+b[i]*n[j]+n[i]*b[j];
    double s=std::clamp((r/H-.4)/.35,0.,1.);d.chi=s*s*s*(10+s*(-15+6*s));d.dchi=(s>0&&s<1)?30*s*s*(1-s)*(1-s)/(.35*H):0;return d;
}
void Count(WallEvaluationStats*stats,const FastWJStatus&status)
{if(stats){stats->fallbacks+=status.UsedFallback()?1:0;stats->rootIterations+=status.rootIterations;stats->maximumRootResidual=std::max(stats->maximumRootResidual,status.normalizedRootResidual);}}
WallStress Value(const Data&d,double k,const Tensor&g,double nu,double a,const SolverUtils::RANSClosureValues&w)
{
    WallStress out;out.chi=d.chi;out.mu=nu+(1-d.chi)*a+d.chi*w.eddyViscosity;
    Require(out.mu>0,"v2 nonpositive implicit diffusion");
    for(int i=0;i<3;++i)for(int j=0;j<3;++j){out.R[i][j]=(1-d.chi)*d.iso[i][j]+d.chi*w.reynoldsStress[i][j];out.V[i][j]=nu*d.strain[i][j]-out.R[i][j]+(i==j?2*k/3:0);out.remainder[i][j]=out.V[i][j]-out.mu*d.strain[i][j];out.production-=out.R[i][j]*g[i][j];}return out;
}
}
WallStress AdaptAmended(double k,double e,const Tensor&g,const Vector&n,double nu,double a,double r,double H,ConstitutiveBackend backend,WallEvaluationStats*stats)
{
    const auto d=Inputs(k,e,g,n,nu,a,r,H);SolverUtils::RANSClosureValues w{};
    if(d.chi>0){if(stats)++stats->values;
        if(backend==ConstitutiveBackend::FastFull3DV1){FastWJStatus status;w=EvaluateFastWJClosure(k,e,g,{},&status);Count(stats,status);}
        else w=SolverUtils::EvaluateWJClosure(k,e,g,{},SolverUtils::WJGradientPolicy::DevStrainRawWorkV1,0);
    }
    return Value(d,k,g,nu,a,w);
}
WallStress ApplyAmendedDirection(double k,double e,const Tensor&g,const Vector&n,double nu,double a,double r,double H,double dk,double de,const Tensor&dg,double da,double dr,ConstitutiveBackend backend,WallEvaluationStats*stats)
{
    Require(std::isfinite(dk)&&std::isfinite(de)&&std::isfinite(da)&&std::isfinite(dr),"nonfinite v2 wall direction");const auto d=Inputs(k,e,g,n,nu,a,r,H);
    Tensor dm{};Vector dmn{},db{};
    for(int i=0;i<3;++i)for(int j=0;j<3;++j){Require(std::isfinite(dg[i][j]),"nonfinite v2 gradient direction");dm[i][j]=dg[i][j]+dg[j][i];dmn[i]+=dm[i][j]*n[j];}
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)db[i]-=d.projector[i][j]*(da*d.mn[j]+a*dmn[j]);
    SolverUtils::RANSClosureValues w{},dw{};
    if(d.chi>0){if(stats){++stats->values;++stats->directions;}
        if(backend==ConstitutiveBackend::FastFull3DV1){auto pair=EvaluateFastWJValueAndDirection(k,e,g,dk,de,dg);w=pair.value;dw=pair.direction;Count(stats,pair.status);}
        else{w=SolverUtils::EvaluateWJClosure(k,e,g,{},SolverUtils::WJGradientPolicy::DevStrainRawWorkV1,0);dw=SolverUtils::ApplyWJClosureDirection(k,e,g,dk,de,dg,{},SolverUtils::WJGradientPolicy::DevStrainRawWorkV1,0);}
    }
    const auto value=Value(d,k,g,nu,a,w);
    const double dc=d.dchi*dr;WallStress out;out.chi=dc;out.mu=(1-d.chi)*da+d.chi*dw.eddyViscosity+dc*(w.eddyViscosity-a);
    for(int i=0;i<3;++i)for(int j=0;j<3;++j){double di=(i==j?2*dk/3:0)+db[i]*n[j]+n[i]*db[j];out.R[i][j]=(1-d.chi)*di+d.chi*dw.reynoldsStress[i][j]+dc*(w.reynoldsStress[i][j]-d.iso[i][j]);out.V[i][j]=nu*dm[i][j]-out.R[i][j]+(i==j?2*dk/3:0);out.remainder[i][j]=out.V[i][j]-out.mu*d.strain[i][j]-value.mu*dm[i][j];out.production-=out.R[i][j]*g[i][j]+value.R[i][j]*dg[i][j];}return out;
}
}
