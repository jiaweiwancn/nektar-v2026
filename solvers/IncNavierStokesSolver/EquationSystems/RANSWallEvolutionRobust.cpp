///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <cmath>
#include <fstream>
#include <cstring>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>

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
    Require(!AcceptedAWF() || !m_cappedDestruction,"accepted budget AWF has no capped-destruction normalization");
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
    std::string reuse;session->LoadSolverInfo("RANSWallFiniteReferenceReuse",reuse,"Off");
    Require(reuse=="Off" || reuse=="ExactHeightV1","unknown immutable reference reuse");
    m_finiteReferenceReuse=reuse=="ExactHeightV1";
    std::string mode;session->LoadSolverInfo("RANSWallRegularizedReference",mode,"NeutralLogV1");
    Require(mode=="NeutralLogV1" || mode=="FiniteInflowV1","unknown prescribed Regularized reference mode");
    if(mode=="FiniteInflowV1")
    {
        Require(AcceptedAWF() && m_variantName=="Regularized" && m_regularizationGamma==0,
            "FiniteInflowV1 requires accepted budget AWF, Regularized CLS and gamma0=0");
        std::string file,expected;session->LoadSolverInfo("RANSWallReferenceFile",file,"");
        session->LoadSolverInfo("RANSWallReferenceSHA256",expected,"");
        Require(!file.empty() && expected.size()==64 && RANSWallR3::HashFile(file)==expected,
            "finite-reference content SHA256 differs or is absent");
        std::ifstream in(file);std::string text((std::istreambuf_iterator<char>(in)),{});
        RANSVector origin{},flow{1,0,0},normal{0,0,1};
        for(int a=0;a<3;++a)
        {
            const std::string axis=a==0?"X":(a==1?"Y":"Z");
            session->LoadParameter("RANSWallReferenceOrigin"+axis,origin[a],0.);
            session->LoadParameter("RANSWallReferenceFlow"+axis,flow[a],a==0?1.:0.);
            session->LoadParameter("RANSWallReferenceNormal"+axis,normal[a],a==2?1.:0.);
        }
        Require(Mapped() || (origin[2]==0. && normal==RANSVector({0,0,1}) && flow==RANSVector({1,0,0})),
            "affine finite initialization/readback requires zero normal origin and horizontal incoming frame");
        m_finiteInflow=std::make_shared<const RANSFiniteInflowReference>(
            RANSFiniteInflowReference::Parse(text),origin,flow,normal,m_regularizationScale,m_parameters.model.Parameters());
        m_referenceFlow=flow;m_referenceNormal=normal;
        m_finiteInflowIdentity=m_finiteInflow->Interpretation()+":content-sha256="+expected;
    }

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
    if(m_assertKCoverage)Require(m_coveredKLocations.count(p),"unprotected constitutive k evaluation: outer closure");
    if(m_variantName=="Original")return EvaluateRANSClosure(m_parameters.model.Model(),k,epsilon,
        gradient,Distance(p),m_parameters.model.Parameters(),0);
    const auto variant=RegularizedReferenceAt(p);
    return EvaluateRANSRobustClosure(m_parameters.model.Model(),k,epsilon,gradient,Distance(p),
        m_parameters.model.Parameters(),0,variant);
}

RANSRobustVariant RANSWallEvolutionPatch::RegularizedReferenceAt(const RANSWallCoordinate &p) const
{
    auto variant=m_robust;
    if(m_finiteInflow)
    {
        const double z=m_finiteInflow->Height(Physical(p));std::uint64_t key;std::memcpy(&key,&z,sizeof key);
        RANSFiniteInflowReference::Value fixed;
        const auto old=m_fixedReferenceValues.find(key);
        if(m_finiteReferenceReuse && old!=m_fixedReferenceValues.end()){fixed=old->second;++m_fixedReferenceHits;}
        else
        {
            fixed=m_finiteInflow->AtHeight(z);++m_fixedReferenceMisses;
            // This caches the immutable mathematical height function, NOT a
            // native field/gradient row. Geometry/side cannot affect it once
            // the physical incoming height is formed. Each patch owns its
            // fixed content/frame/interpretation; exact height bits distinguish
            // signed zero. Moving heights cannot grow storage without bound.
            if(m_finiteReferenceReuse)
            {if(m_fixedReferenceValues.size()>=32768){m_fixedReferenceValues.clear();++m_fixedReferenceClears;}
             m_fixedReferenceValues.emplace(key,fixed);}
        }
        variant.gamma=fixed.gamma;variant.referenceStrain=fixed.strain;return variant;
    }
    // Preserve the historical arithmetic and default identity exactly.
    const double r=Distance(p),kap=m_parameters.model.Kappa();
    variant.gamma=m_regularizationGamma+m_regularizationScale*kap*m_robustUstar*r;
    const double shear=m_robustUstar/(kap*r);
    for(int a=0;a<3;++a)for(int b=0;b<3;++b)
        variant.referenceStrain[a][b]=shear*(m_referenceFlow[a]*m_referenceNormal[b]+m_referenceNormal[a]*m_referenceFlow[b]);
    return variant;
}

RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::FiniteInflowSeed() const
{
    Require(m_finiteInflow && !Mapped() && AcceptedAWF(),"finite-profile seed currently requires the affine accepted-AWF pilot");
    const double H=m_parameters.geometry.height;
    const auto endpoint=m_finiteInflow->AtHeight(H);const auto wall=Profile(endpoint.k,{0,0,0},1);
    Require(m_dimension==2 && m_referenceFlow==RANSVector({1,0,0}) && m_referenceNormal==RANSVector({0,0,1}),
        "finite compatibility seed requires the horizontal Cartesian incoming frame");
    std::vector<double> breaks{wall->Reference().cutoff,wall->Original().ViscosityInterceptHeight(),m_parameters.geometry.roughness};
    for(const auto &q:m_finiteInflow->Knots())if(q.z>0 && q.z<upper[2])breaks.push_back(q.z);
    const auto velocity=[&](double r){return r<H?endpoint.U*wall->Original().Integrals(r).j0/wall->Original().Integrals().j0:m_finiteInflow->AtHeight(r).U;};
    Vector s(Size());
    const auto u=fields[0]->RadialProject(velocity,{{0,0,0},{H,0,endpoint.U},{upper[2],0,velocity(upper[2])}},breaks);
    std::copy(u.begin(),u.end(),s.begin()+offset[0]);
    const auto k=fields[K()]->RadialProject([&](double r){return r<H?wall->ReferenceK(r):m_finiteInflow->AtHeight(r).k;},
        {{0,1,0},{H,0,endpoint.k},{upper[2],0,m_finiteInflow->AtHeight(upper[2]).k}},breaks);
    std::copy(k.begin(),k.end(),s.begin()+offset[K()]);
    const auto epsilon=fields[Epsilon()]->RadialProject([&](double r){return m_finiteInflow->AtHeight(r).epsilon;},
        {{H,0,wall->Dissipation(H)},{upper[2],0,m_finiteInflow->AtHeight(upper[2]).epsilon}},breaks);
    std::copy(epsilon.begin(),epsilon.end(),s.begin()+offset[Epsilon()]);
    return s;
}

std::string RANSWallEvolutionPatch::FiniteInflowChecks(const Vector &s) const
{
    namespace J=RANSWallSolveJson;
    Require(bool(m_finiteInflow),"finite-reference diagnostic needs its explicit provider");
    double workError=0,epsilonError=0,trace=0,symmetry=0,minGamma=1e100,maxGamma=0;
    double minK=1e100,minE=1e100;int points=0;
    for(const auto &q:m_outer)
    {
        const auto fixed=m_finiteInflow->At(Physical(q.v.coordinate));
        minGamma=std::min(minGamma,fixed.gamma);maxGamma=std::max(maxGamma,fixed.gamma);
        const auto v=RegularizedReferenceAt(q.v.coordinate);
        Require(v.gamma==fixed.gamma && v.referenceStrain==fixed.strain,"provider/assembly reference differs");
        trace=std::max(trace,std::abs(v.referenceStrain[0][0]+v.referenceStrain[1][1]+v.referenceStrain[2][2]));
        for(int a=0;a<3;++a)for(int b=0;b<3;++b)symmetry=std::max(symmetry,std::abs(v.referenceStrain[a][b]-v.referenceStrain[b][a]));
        const double k=q.v.Value(s,offset[K()]),e=q.e.Value(s,offset[Epsilon()]);
        minK=std::min(minK,k);minE=std::min(minE,e);
        if(k>0 && e>0)
        {
            RANSTensor g{};for(int a=0;a<m_dimension;++a)g[PhysicalAxis(a)]=q.v.Gradient(s,offset[a]);
            const auto c=EvaluateRANSRobustClosure(m_parameters.model.Model(),k,e,g,Distance(q.v.coordinate),m_parameters.model.Parameters(),0,v);
            long double work=0;for(int a=0;a<3;++a)for(int b=0;b<3;++b)work-=(long double)c.reynoldsStress[a][b]*g[a][b];
            workError=std::max(workError,std::abs(double(work)-c.production)/std::max(1.,std::abs(c.production)));
            epsilonError=std::max(epsilonError,std::abs(c.epsilonProduction-m_parameters.model.Parameters().cEpsilon1*e/k*c.production)/std::max(1.,std::abs(c.epsilonProduction)));
        }
        ++points;
    }
    std::vector<std::string> samples;
    for(double z:{m_parameters.geometry.height,.04,.1,.25,.4,.8,1.})
    {
        const auto v=m_finiteInflow->AtHeight(z);
        samples.push_back(J::Numbers(Vector{z,v.U,v.dU,v.k,v.epsilon,v.nuT,v.gamma}));
    }
    return J::Object({{"reference_reuse_mode",J::Quote(m_finiteReferenceReuse?"ExactHeightV1":"Off")},
        {"reference_reuse_hits_misses_entries_clears",J::Numbers(Vector{double(m_fixedReferenceHits),double(m_fixedReferenceMisses),double(m_fixedReferenceValues.size()),double(m_fixedReferenceClears)})},
        {"mode",J::Quote("FiniteInflowV1")},{"identity",J::Quote(m_finiteInflowIdentity)},
        {"comparison_reference",J::Quote("native finite-profile initialized seed; not analytical neutral equilibrium")},
        {"frame_height",J::Quote("dot(physical-X_in,n_in); independent of wall distance and outer chart")},
        {"outer_native_points",J::Number(points)},{"gamma_min_max",J::Numbers(Vector{minGamma,maxGamma})},
        {"strain_trace_symmetry_defects",J::Numbers(Vector{trace,symmetry})},
        {"stress_work_relative_defect",J::Number(workError)},{"epsilon_production_relative_defect",J::Number(epsilonError)},
        {"outer_native_scalar_minima",J::Numbers(Vector{minK,minE})},
        {"sample_columns",J::Quote("zeta,U_ref,U_ref_prime,k_ref,epsilon_ref,nu_t_ref,gamma")},
        {"samples",J::List(samples)}});
}

RANSCappedDestructionValues RANSWallEvolutionPatch::WallNormalization(const Vector &s,
    const RANSWallCoordinate &p,const EquilibriumAWFProfiles &profile) const
{
    Require(!AcceptedAWF(),"capped normalization is not part of CLSCalibratedBudgetAWF");
    if(auto cached=CachedNormalization(s,p,profile))return *cached;
    ++m_evaluationWork.normalizationComputes;
    RANSCappedDestructionInput input;input.sampledK=profile.K();
    // Every rank evaluates the complete ray. Ownership is applied only when
    // assembling its volume loads, never to the nonlocal normalization.
    for(const auto &q:profile.Quadrature())
    {
        auto point=p;point[2]=q.r/m_mapping.At(p).scale;
        const auto row=EvaluationRow(K(),fields[K()]->LocateElement(point),point);
        if(m_assertKCoverage)Require(m_coveredKLocations.count(point),
            "unprotected constitutive k evaluation: normalization ray");
        const double k=row->Value(s,offset[K()]);input.k.push_back(k);
        if(m_solveInvalidContext.empty() && (!std::isfinite(k) || k<0))
            m_solveInvalidContext=RANSWallSolveJson::Object({
                {"operation",RANSWallSolveJson::Quote("complete_ray_capped_normalization")},
                {"chart",RANSWallSolveJson::Numbers(point)},{"local_k",RANSWallSolveJson::Number(k)},
                {"geometry",RANSWallSolveJson::Number(row->geometry)},
                {"sampled_K",RANSWallSolveJson::Number(profile.K())},
                {"ray_chart",RANSWallSolveJson::Numbers(p)}});
        Require(std::isfinite(k) && k>=0,"Capped destruction requires finite nonnegative transported k");
    }
    for(const auto &q:profile.Quadrature())
    {
        const auto pre=profile.AtProduction(q.r,0);
        input.preE.push_back(pre.dissipation);
        input.preD.push_back(pre.dissipation-pre.redistribution);
        input.metricWeights.push_back(q.weight*q.jacobian);
    }
    auto result=EvaluateRANSCappedDestruction(input,m_destruction);
    StoreNormalization(s,p,profile,result);return result;
}

AWFProfileValues RANSWallEvolutionPatch::WallSource(const Vector &s,const RANSWallCoordinate &p,
    double production,const EquilibriumAWFProfiles &profile,const RANSCappedDestructionValues *normalization,
    double *destructionRate) const
{
    auto values=profile.AtProduction(Distance(p),production);
    if(!m_cappedDestruction)return values;
    RANSCappedDestructionValues local;
    if(!normalization){local=WallNormalization(s,p,profile);normalization=&local;}
    if(m_assertKCoverage)Require(m_coveredKLocations.count(p),"unprotected constitutive k evaluation: wall source");
    const double k=Point(K(),p).Value(s,offset[K()]);
    const auto repaired=EvaluateRANSCappedDestructionPoint(k,profile.K(),values.dissipation,
        values.dissipation-values.redistribution,normalization->b,normalization->m,m_destruction);
    if(destructionRate)*destructionRate=repaired.rate;
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
