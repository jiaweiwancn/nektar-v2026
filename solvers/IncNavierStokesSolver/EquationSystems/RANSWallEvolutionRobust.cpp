///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <cmath>

namespace Nektar
{
using namespace SolverUtils;
using RANSWallAlgebra::Require;

void RANSWallEvolutionPatch::InitialiseRobust(const LibUtilities::SessionReaderSharedPtr &session)
{
    std::string innerStress;session->LoadSolverInfo("RANSWallInnerStress",innerStress,"Linear");
    m_innerStress=ParseRANSWallInnerStress(innerStress);
    session->LoadSolverInfo("RANSWallVariant",m_variantName,"Original");
    Require(m_variantName=="Original" || m_variantName=="CubicScaled" || m_variantName=="Regularized",
        "RANSWallVariant must be Original, CubicScaled or Regularized");
    const bool selected=m_variantName!="Original";
    session->LoadParameter("RANSWallCubicBeta",m_robust.beta,selected?.4:1.);
    session->LoadParameter("RANSWallRegularizationGamma",m_regularizationGamma,0.);
    session->LoadParameter("RANSWallRegularizationScale",m_regularizationScale,0.);
    Require(std::isfinite(m_regularizationGamma) && m_regularizationGamma>=0 &&
        std::isfinite(m_regularizationScale) && m_regularizationScale>=0,"invalid prescribed regularization viscosity");
    Require(m_variantName=="Regularized" || (m_regularizationGamma==0 && m_regularizationScale==0),
        "positive increment requires the explicit Regularized model selection");
    Require(selected || m_robust.beta==1,"Original selects the frozen unmodified cubic coefficients");
    Require(!selected || m_parameters.model.Model()==RANSModel::CLSHighRe,"CLS variants require CLSHighRe");
    m_robust.Validate();
    std::string destruction,smoothing;
    session->LoadSolverInfo("RANSWallDestruction",destruction,selected?"Capped":"Original");
    Require(destruction=="Original" || destruction=="Capped","unknown inner destruction selection");
    m_cappedDestruction=destruction=="Capped";
    session->LoadSolverInfo("RANSWallDestructionSmoothing",smoothing,"Linear");
    Require(smoothing=="Linear" || smoothing=="C2","unknown destruction limiter shape");
    m_destruction.smooth=smoothing=="C2";
    session->LoadParameter("RANSWallPositivityTheta",m_destruction.theta,.2);
    session->LoadParameter("RANSWallPositivityDelta",m_destruction.delta,.1);
    m_destruction.Validate();
    session->LoadParameter("RANSWallFrictionVelocity",m_robustUstar,.4);
    double dx,dy;session->LoadParameter("RANSWallDirectionX",dx,1.);session->LoadParameter("RANSWallDirectionY",dy,0.);
    const double length=std::hypot(dx,dy);Require(length>0 && m_robustUstar>0,"invalid regularization reference inputs");
    const auto incoming=m_mapping.At({lower[0],lower[1],0});
    m_referenceNormal=incoming.normal;
    for(int a=0;a<3;++a)m_referenceFlow[a]=(dx*incoming.tangent1[a]+dy*incoming.tangent2[a])/length;
}

std::string RANSWallEvolutionPatch::InnerStressName() const
{
    return RANSWallInnerStressName(m_innerStress);
}

RANSWallInnerStressValues RANSWallEvolutionPatch::InnerStress(double k,double auxiliaryNu,
    const RANSTensor &gradient) const
{
    return EvaluateRANSWallInnerStress(k,auxiliaryNu,gradient,m_innerStress);
}

RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::RobustParameters() const
{
    return {m_robust.beta,m_regularizationGamma,m_destruction.theta,m_destruction.delta,
        m_destruction.smooth?1.:0.,m_regularizationScale,m_cappedDestruction?1.:0.};
}

RANSClosureValues RANSWallEvolutionPatch::OuterClosure(double k,double epsilon,
    const RANSTensor &gradient,const RANSWallCoordinate &p) const
{
    if(m_variantName=="Original")return EvaluateRANSClosure(m_parameters.model.Model(),k,epsilon,
        gradient,Distance(p),m_parameters.model.Parameters(),0);
    auto variant=m_robust;
    const double r=Distance(p),kap=m_parameters.model.Kappa();
    variant.gamma=m_regularizationGamma+m_regularizationScale*kap*m_robustUstar*r;
    // Prescribed Cartesian extension of the incoming flat strain. It uses
    // chart-normal physical distance and the frozen incoming frame, with no
    // dependence on the evolving velocity, k, epsilon, or hill normal.
    const double shear=m_robustUstar/(kap*r);
    for(int a=0;a<3;++a)for(int b=0;b<3;++b)
        variant.referenceStrain[a][b]=shear*(m_referenceFlow[a]*m_referenceNormal[b]+
            m_referenceNormal[a]*m_referenceFlow[b]);
    return EvaluateRANSRobustClosure(m_parameters.model.Model(),k,epsilon,gradient,r,
        m_parameters.model.Parameters(),0,variant);
}

RANSCappedDestructionValues RANSWallEvolutionPatch::WallNormalization(const Vector &s,
    const RANSWallCoordinate &p,const EquilibriumAWFProfiles &profile) const
{
    RANSCappedDestructionInput input;input.sampledK=profile.K();
    const double scale=m_mapping.At(p).scale;
    // Every rank evaluates the complete ray. Ownership is applied only when
    // assembling its volume loads, never to the nonlocal normalization.
    for(const auto &q:profile.Quadrature())
    {
        auto point=p;point[2]=q.r/scale;
        const auto pre=profile.AtProduction(q.r,0);
        const double k=Point(K(),point).Value(s,offset[K()]);input.k.push_back(k);
        if(m_solveDiagnostics && m_solveInvalidContext.empty() && (!std::isfinite(k) || k<0))
            m_solveInvalidContext=RANSWallSolveJson::Object({
                {"operation",RANSWallSolveJson::Quote("complete_ray_capped_normalization")},
                {"chart",RANSWallSolveJson::Numbers(point)},{"local_k",RANSWallSolveJson::Number(k)},
                {"sampled_K",RANSWallSolveJson::Number(profile.K())},
                {"ray_chart",RANSWallSolveJson::Numbers(p)}});
        input.preE.push_back(pre.dissipation);
        input.preD.push_back(pre.dissipation-pre.redistribution);
        input.metricWeights.push_back(q.weight*q.jacobian);
    }
    return EvaluateRANSCappedDestruction(input,m_destruction);
}

AWFProfileValues RANSWallEvolutionPatch::WallSource(const Vector &s,const RANSWallCoordinate &p,
    double production,const EquilibriumAWFProfiles &profile,const RANSCappedDestructionValues *normalization) const
{
    auto values=profile.AtProduction(Distance(p),production);
    if(!m_cappedDestruction)return values;
    RANSCappedDestructionValues local;
    if(!normalization){local=WallNormalization(s,p,profile);normalization=&local;}
    const double k=Point(K(),p).Value(s,offset[K()]);
    const auto repaired=EvaluateRANSCappedDestructionPoint(k,profile.K(),values.dissipation,
        values.dissipation-values.redistribution,normalization->b,normalization->m,m_destruction);
    values.dissipation=repaired.dissipation;values.redistribution=repaired.redistribution;
    values.kSource=production-repaired.destruction;
    return values;
}

double RANSWallEvolutionPatch::WallDissipation(const Vector &s,const RANSWallCoordinate &p,
    const EquilibriumAWFProfiles &profile) const
{
    if(!m_cappedDestruction)return profile.Dissipation(Distance(p));
    return WallNormalization(s,p,profile).lambda*profile.Dissipation(Distance(p));
}
}
