///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <SolverUtils/Turbulence/RANSWallInnerStress.h>
#include <cmath>
#include <stdexcept>

namespace Nektar::SolverUtils
{
namespace
{
void Check(bool valid,const char *message)
{
    if(!valid)throw std::invalid_argument(message);
}
void Finite(const RANSTensor &a)
{
    for(const auto &row:a)for(double value:row)
        Check(std::isfinite(value),"nonfinite inner-stress tensor input");
}
double Norm(const RANSTensor &a)
{
    double answer=0;for(const auto &row:a)for(double value:row)answer=std::hypot(answer,value);
    return answer;
}
}

const char *RANSWallInnerStressName(RANSWallInnerStressModel model)
{
    switch(model)
    {
        case RANSWallInnerStressModel::Linear:return "Linear";
        case RANSWallInnerStressModel::SaturatedR2:return "SaturatedR2";
    }
    throw std::invalid_argument("unknown inner-stress model enum");
}
RANSWallInnerStressModel ParseRANSWallInnerStress(const std::string &name)
{
    if(name=="Linear")return RANSWallInnerStressModel::Linear;
    if(name=="SaturatedR2")return RANSWallInnerStressModel::SaturatedR2;
    throw std::invalid_argument("RANSWallInnerStress must be Linear or SaturatedR2");
}

RANSWallInnerStressValues EvaluateRANSWallInnerStress(double k,double auxiliaryNu,
    const RANSTensor &gradient,RANSWallInnerStressModel model)
{
    Check(std::isfinite(k) && k>=0,"negative/nonfinite local k in inner stress");
    Check(std::isfinite(auxiliaryNu) && auxiliaryNu>=0,"negative/nonfinite auxiliary viscosity in inner stress");
    RANSWallInnerStressName(model);Finite(gradient);
    RANSWallInnerStressValues value;value.model=model;value.k=k;value.auxiliaryNu=auxiliaryNu;
    value.c=2*k/3;
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)value.strain[i][j]=gradient[i][j]+gradient[j][i];
    Finite(value.strain);value.strainNorm=Norm(value.strain);value.a=auxiliaryNu*value.strainNorm;
    Check(std::isfinite(value.c) && std::isfinite(value.a),"unrepresentable inner-stress scale");
    value.ta=1;value.radialTangent=value.tangentialTangent=auxiliaryNu;value.t=value.a;
    if(model==RANSWallInnerStressModel::SaturatedR2 && k==0)
    {
        value.linearResponse=false;value.t=0;value.ta=0;
        value.radialTangent=value.tangentialTangent=0;
    }
    else if(model==RANSWallInnerStressModel::SaturatedR2 && value.a>RANSWallInnerSaturationRho*value.c)
    {
        const double d=1-RANSWallInnerSaturationRho;
        const double z=(value.a-RANSWallInnerSaturationRho*value.c)/(d*value.c);
        const double q=std::exp(-z);
        value.active=true;value.linearResponse=false;value.exponentialUnderflow=q==0;
        value.ta=q;
        value.tc=q==0?1:-std::expm1(-z)-d*z*q;
        value.t=value.c*(RANSWallInnerSaturationRho-d*std::expm1(-z));
        value.radialTangent=auxiliaryNu*q;
        value.tangentialTangent=value.t/value.strainNorm;
    }
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)
    {
        value.turbulentStress[i][j]=value.linearResponse?auxiliaryNu*value.strain[i][j]:
            (value.strainNorm>0?value.t*(value.strain[i][j]/value.strainNorm):0);
        value.reynoldsStress[i][j]=(i==j?value.c:0)-value.turbulentStress[i][j];
        value.production-=value.reynoldsStress[i][j]*gradient[i][j];
    }
    Check(std::isfinite(value.production),"nonfinite inner-stress production");
    return value;
}

RANSWallInnerStressAction RANSWallInnerStressDirectional(double k,double auxiliaryNu,
    const RANSTensor &gradient,double dk,double dn,const RANSTensor &h,RANSWallInnerStressModel model)
{
    const auto value=EvaluateRANSWallInnerStress(k,auxiliaryNu,gradient,model);
    Check(std::isfinite(dk) && std::isfinite(dn),"nonfinite inner-stress scalar direction");Finite(h);
    Check(k>0 || dk>=0,"inner-stress direction leaves nonnegative local-k cone");
    Check(auxiliaryNu>0 || dn>=0,"inner-stress direction leaves nonnegative viscosity cone");
    RANSWallInnerStressAction action;
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)action.strain[i][j]=h[i][j]+h[j][i];
    Finite(action.strain);
    if(model==RANSWallInnerStressModel::Linear || value.linearResponse)
    {
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)
            action.turbulentStress[i][j]=auxiliaryNu*action.strain[i][j]+dn*value.strain[i][j];
    }
    else if(k>0)
    {
        double projection=0;
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)
            projection+=(value.strain[i][j]/value.strainNorm)*action.strain[i][j];
        const double radial=(value.radialTangent-value.tangentialTangent)*projection;
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)
            action.turbulentStress[i][j]=value.tangentialTangent*action.strain[i][j]+
                radial*(value.strain[i][j]/value.strainNorm)+value.ta*dn*value.strain[i][j]+
                value.tc*(2*dk/3)*(value.strain[i][j]/value.strainNorm);
    }
    else if(auxiliaryNu>0 && value.strainNorm>0)
    {
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)
            action.turbulentStress[i][j]=(2*dk/3)*(value.strain[i][j]/value.strainNorm);
    }
    else if(auxiliaryNu>0 && value.strainNorm==0)
    {
        action.directionalOnly=true;
        action.turbulentStress=EvaluateRANSWallInnerStress(dk,auxiliaryNu,h,model).turbulentStress;
    }
    else if(auxiliaryNu==0 && value.strainNorm>0)
    {
        action.directionalOnly=true;
        action.turbulentStress=EvaluateRANSWallInnerStress(dk,dn,gradient,model).turbulentStress;
    }
    // At k=nuA=S=0, F is second order and its Frechet action is zero.
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)
    {
        action.reynoldsStress[i][j]=(i==j?2*dk/3:0)-action.turbulentStress[i][j];
        action.production-=action.reynoldsStress[i][j]*gradient[i][j]+value.reynoldsStress[i][j]*h[i][j];
    }
    Finite(action.turbulentStress);Finite(action.reynoldsStress);
    Check(std::isfinite(action.production),"nonfinite inner-stress production action");
    return action;
}

RANSTensor RANSWallInnerMeanStress(const RANSWallInnerStressValues &value,double nu,double pi)
{
    Check(std::isfinite(nu) && nu>=0 && std::isfinite(pi),"invalid inner mean-stress inputs");
    RANSTensor stress{};
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)
        stress[i][j]=(value.linearResponse?(nu+value.auxiliaryNu)*value.strain[i][j]:
            nu*value.strain[i][j]+value.turbulentStress[i][j])-(i==j?pi:0);
    Finite(stress);return stress;
}
} // namespace Nektar::SolverUtils
