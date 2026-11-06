///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <algorithm>
#include <cmath>

namespace Nektar
{
using namespace SolverUtils;
using RANSWallAlgebra::Require;
namespace
{
bool Same(double a,double b){return std::abs(a-b)<1e-11;}
double Dot(const RANSVector &a,const RANSVector &b)
{return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
}

double RANSWallEvolutionPatch::WallCoordinate(const RANSWallCoordinate &p) const
{return EquilibriumAWFWallCoordinate(m_mapping,p);}
double RANSWallEvolutionPatch::DependentEpsilon(const Vector &s,const RANSWallCoordinate &p,
    double time,const RANSWallEvolutionData &data) const
{return WallDissipation(s,p,*StateProfile(s,p,time,data));}

void RANSWallEvolutionPatch::InitialiseTransition(const LibUtilities::SessionReaderSharedPtr &session)
{
    Require(m_parameters.chi==0 || m_parameters.chi==1,"Upstream selects the prescribed weight (chi=1) or disables it (chi=0)");
    session->LoadParameter("RANSWallTransitionStart",m_transition.start,-3.);
    session->LoadParameter("RANSWallTransitionEnd",m_transition.end,-2.25);
    session->LoadParameter("RANSWallBaselineStart",m_baselineStart,-2.);
    session->LoadParameter("RANSWallSensorDistanceA",m_sensorA,2.5);
    session->LoadParameter("RANSWallSensorDistanceB",m_sensorB,3.5);
    auto &c=m_transition.state;
    for(auto entry:std::vector<std::pair<std::string,double*>>{{"On",&c.on},{"Off",&c.off},
        {"MarginOff",&c.marginOff},{"MarginOn",&c.marginOn},{"ReynoldsOff",&c.reynoldsOff},{"ReynoldsOn",&c.reynoldsOn}})
        session->LoadParameter("RANSWallSensor"+entry.first,*entry.second,*entry.second);
    const std::array<std::string,7> names{"Velocity","Increment","Normal","K","Epsilon","Pressure","Convection"};
    for(int j=0;j<7;++j)session->LoadParameter("RANSWallSensorWeight"+names[j],c.weights[j],1.);
    m_transition.Validate();
    Require(c.reynoldsOff>=m_options.minimumReynolds,"the smooth Reynolds gate must stay inside the declared high-Re branch");
    if(!m_uniformEP)ValidateEquilibriumAWFTransitionSupport(m_mapping,m_transition);
    Require(m_mapping.heightSlopeX==0 && m_mapping.heightSlopeY==0 && m_mapping.heightAmplitude==0,
        "the flat upstream reference requires constant physical H; use the adapted baseline otherwise");
    Require(m_uniformEP || (std::isfinite(m_baselineStart) && m_baselineStart>m_transition.end),
        "the designated baseline faces need a positive physical transition buffer");
    if(m_mapping.kind==RANSWallNormalMap::Kind::CosineHill)
        Require(m_baselineStart<WallCoordinate({-m_mapping.lengthX,0,0}),"the projection buffer must end before the incoming hill foot");
    Require(m_sensorA>m_mapping.height && m_sensorA<m_sensorB && m_sensorB<upper[2],
        "independent sensors require H<r_a<r_b<outer top");
    const auto atH=[&](int axis,int side,const RANSWallGlobalField::Element &e)
        {return axis==2 && !side && Same(e.lower[2],m_parameters.geometry.height);};
    for(auto [id,value]:fields[Epsilon()]->Boundary([](const RANSWallCoordinate&){return 0.;},atH))
    {(void)value;m_traceIds.insert(id);}
    if(m_uniformEP)
    {
        // No downstream baseline-only trace modes in the uniform control.
        // Native reference/essential data and initialization use this same
        // fully active profile; no old state or history is migrated.
        Require(m_dimension==2 && !Mapped() && m_parameters.chi==1,"UniformEPFlat scope");
        double us,dx,dy;session->LoadParameter("RANSWallFrictionVelocity",us,.4);
        session->LoadParameter("RANSWallDirectionX",dx,1.);session->LoadParameter("RANSWallDirectionY",dy,0.);
        m_transitionReference=FlatReference(us,{dx,dy});m_transitionReady=true;
        const auto profile=Profile(m_parameters.model.Q()*us*us,{lower[0],0,0},1);
        Require(profile->Reference().weight==1 && profile->Reference().reason==AWFFallbackReason::Active,
            "UniformEPFlat reference is not intrinsically admissible");
        return;
    }
    bool plateau=false,baseline=false;
    for(const auto &e:fields[Epsilon()]->Elements())if(Same(e.lower[2],m_parameters.geometry.height))
    {
        const double a=WallCoordinate(e.lower),b=WallCoordinate(e.upper);
        Require(!(a<m_transition.start-1e-11 && b>m_transition.start+1e-11),
            "a native matching face must end at transition start to preserve the flat plateau");
        Require(!(a<m_baselineStart-1e-11 && b>m_baselineStart+1e-11),
            "designated baseline support must start on a native matching-face boundary");
        plateau=plateau || b<=m_transition.start+1e-11;
        baseline=baseline || a>=m_baselineStart-1e-11;
    }
    Require(plateau && baseline,"Upstream needs both a flat plateau and designated baseline faces");
    const auto baselineFace=[&](int axis,int side,const RANSWallGlobalField::Element &e)
        {return atH(axis,side,e) && WallCoordinate(e.lower)>=m_baselineStart-1e-11;};
    for(auto [id,value]:fields[Epsilon()]->Boundary([](const RANSWallCoordinate&){return 0.;},baselineFace))
    {(void)value;m_baselineTraceIds.insert(id);}
    double us,dx,dy;session->LoadParameter("RANSWallFrictionVelocity",us,.4);
    session->LoadParameter("RANSWallDirectionX",dx,1.);session->LoadParameter("RANSWallDirectionY",dy,0.);
    const auto incoming=m_mapping.At({lower[0],lower[1],0});
    const double incomingX=dx*incoming.tangent1[0]+dy*incoming.tangent2[0];
    Require(incomingX>0,"the physical incoming flow must increase the prescribed wall coordinate");
    if(m_mapping.kind==RANSWallNormalMap::Kind::CosineHill)
        Require(dy==0,"the extruded hill transition reference follows its prescribed wall curve");
    m_transitionReference=FlatReference(us,{dx,dy});
    const auto ref=Profile(m_parameters.model.Q()*us*us,{lower[0],lower[1],0},0)->Reference();
    Require(ref.margin>=c.marginOn && ref.reynolds>=c.reynoldsOn,
        "the chosen flat reference must lie on both admissibility unit plateaus");
    m_transitionReady=true;
    // Test the actual stored projection and actual evaluation/forcing maps.
    RANSWallEvolutionData zero(m_dimension);
    for(const auto &ray:m_rays)if(WallCoordinate(ray.coordinate)<m_transition.end)
    {
        const auto result=Transition(m_transitionReference,ray.coordinate,0,zero);
        Require(result.d2<1e-24 && result.state==1,"projected transition reference does not zero every sensor");
    }
}

std::array<RANSVector,2> RANSWallEvolutionPatch::ForcingParts(const Vector &s,
    const RANSWallCoordinate &point,double time,const RANSWallEvolutionData &data) const
{
    auto p=point;p[2]=m_parameters.geometry.height;const auto metric=m_mapping.At(p);
    const auto v=Point(0,p,1),pressure=Point(OuterPressure(),p,1),k=Point(K(),p,1);
    RANSTensor G{},P{},gradient{};RANSVector u{},ut{};
    for(int a=0;a<m_dimension;++a){const int aa=PhysicalAxis(a);u[aa]=v.Value(s,offset[a]);G[aa]=v.Gradient(s,offset[a]);}
    const auto &n=metric.normal;const double w=Dot(u,n);
    for(int a=0;a<3;++a){ut[a]=u[a]-w*n[a];for(int b=0;b<3;++b)P[a][b]=(a==b?1.:0.)-n[a]*n[b];}
    for(int a=0;a<3;++a)for(int b=0;b<3;++b)
    {
        for(int d=0;d<3;++d)gradient[a][b]+=P[a][d]*G[d][b]-n[a]*metric.Br[d][b]*u[d];
        gradient[a][b]-=w*metric.Br[a][b];
    }
    double divergence=0;for(int a=0;a<3;++a)divergence+=gradient[a][a];
    const auto gp=pressure.Gradient(s,offset[OuterPressure()]),gk=k.Gradient(s,offset[K()]);
    RANSVector pressureBody{},convection{};
    for(int a=0;a<3;++a){pressureBody[a]=gp[a]-2*gk[a]/3;convection[a]=ut[a]*divergence;
        for(int b=0;b<3;++b)convection[a]+=gradient[a][b]*ut[b];}
    const auto &body=data.sensorBody.empty()?data.outerBody:data.sensorBody;
    for(int a=0;a<m_dimension;++a)pressureBody[PhysicalAxis(a)]-=body[a](p,time);
    std::array<RANSVector,2> result{};
    for(int a=0;a<3;++a)for(int b=0;b<3;++b)
    {result[0][a]+=P[a][b]*pressureBody[b];result[1][a]+=P[a][b]*convection[b];}
    return result;
}

AWFTransitionInput RANSWallEvolutionPatch::TransitionInput(const Vector &s,
    const RANSWallCoordinate &point,double time,const RANSWallEvolutionData &data) const
{
    const double K0=SampleK(s,point);const auto base=Profile(K0,point,0);
    AWFTransitionInput input;input.coordinate=WallCoordinate(point);
    input.budgetMargin=base->Reference().margin;input.reynolds=base->Reference().reynolds;
    if(input.coordinate>=m_transition.end)return input;
    Require(m_transitionReady,"transition reference has not been prepared");
    auto H=point;H[2]=m_parameters.geometry.height;const auto metric=m_mapping.At(H);
    Require(std::abs(metric.kappa1)+std::abs(metric.kappa2)<1e-14,
        "nonzero spatial correction support must be exactly flat");
    auto a=H,b=H;a[2]=m_sensorA/metric.scale;b[2]=m_sensorB/metric.scale;
    const auto vh=Point(0,H,1),va=Point(0,a),vb=Point(0,b),ka=Point(K(),a),ea=Point(Epsilon(),a);
    AWFTransitionPhysicalInput physical;physical.normal=metric.normal;physical.sampledK=K0;
    physical.height=metric.height;physical.molecularNu=m_parameters.geometry.molecularNu;physical.Q=m_parameters.model.Q();
    for(int j=0;j<m_dimension;++j)
    {
        const int d=PhysicalAxis(j);
        physical.matchingVelocity[d]=vh.Value(s,offset[j]);physical.matchingReference[d]=vh.Value(m_transitionReference,offset[j]);
        physical.velocityA[d]=va.Value(s,offset[j]);physical.velocityAReference[d]=va.Value(m_transitionReference,offset[j]);
        physical.velocityB[d]=vb.Value(s,offset[j]);physical.velocityBReference[d]=vb.Value(m_transitionReference,offset[j]);
    }
    physical.kA=ka.Value(s,offset[K()]);physical.epsilonA=ea.Value(s,offset[Epsilon()]);
    physical.epsilonAReference=ea.Value(m_transitionReference,offset[Epsilon()]);
    // The stored flat family has zero prescribed body force. Evaluate both
    // states with the same discrete forcing operator, retaining that fixed
    // reference input rather than cancelling a new physical body force.
    RANSWallEvolutionData referenceData(m_dimension);
    const auto force=ForcingParts(s,H,time,data),reference=ForcingParts(m_transitionReference,H,time,referenceData);
    physical.pressure=force[0];physical.pressureReference=reference[0];
    physical.convection=force[1];physical.convectionReference=reference[1];
    input.sensors=EvaluateEquilibriumAWFTransitionSensors(physical).value;
    return input;
}

AWFTransitionValues RANSWallEvolutionPatch::Transition(const Vector &s,const RANSWallCoordinate &p,
    double time,const RANSWallEvolutionData &data) const
{
    if(m_uniformEP)
    {
        const auto profile=Profile(SampleK(s,p),p,1);
        Require(profile->Reference().weight==1 && profile->Reference().reason==AWFFallbackReason::Active,
            "UniformEPFlat profile is not intrinsically admissible");
        AWFTransitionValues value;value.chi=value.spatial=value.state=value.admissibility=value.marginGate=value.reynoldsGate=1;
        return value;
    }
    auto result=EvaluateEquilibriumAWFTransition(TransitionInput(s,p,time,data),m_transition).value;
    if(m_parameters.chi==0)result.chi=0;
    return result;
}

std::shared_ptr<EquilibriumAWFProfiles> RANSWallEvolutionPatch::StateProfile(const Vector &s,
    const RANSWallCoordinate &p,double time,const RANSWallEvolutionData &data) const
{
    const double K0=SampleK(s,p);
    if(m_uniformEP)
    {
        const auto profile=Profile(K0,p,1);
        Require(profile->Reference().weight==1 && profile->Reference().reason==AWFFallbackReason::Active,
            "UniformEPFlat forbids an inactive/fallback profile");
        return profile;
    }
    if(!m_transitionEnabled)return Profile(K0,p);
    if(m_parameters.chi==0 || WallCoordinate(p)>=m_transition.end)return Profile(K0,p,0);
    const auto base=Profile(K0,p,0);const auto &ref=base->Reference();
    if(ref.margin<=m_transition.state.marginOff || ref.reynolds<=m_transition.state.reynoldsOff)return base;
    const double weight=Transition(s,p,time,data).chi;
    return weight==0?base:Profile(K0,p,weight);
}

std::array<std::map<int,double>,2> RANSWallEvolutionPatch::DissipationTrace(const Vector &s,
    double time,const RANSWallEvolutionData &data) const
{
    const auto face=[&](int axis,int side,const RANSWallGlobalField::Element &e)
        {return axis==2 && !side && Same(e.lower[2],m_parameters.geometry.height);};
    std::array<std::map<int,double>,2> trace;
    trace[0]=fields[Epsilon()]->Boundary([&](const RANSWallCoordinate &p)
        {return WallDissipation(s,p,*Profile(SampleK(s,p),p,0));},face);
    trace[1]=fields[Epsilon()]->Boundary([&](const RANSWallCoordinate &p)
    {
        const auto baseline=Profile(SampleK(s,p),p,0),current=StateProfile(s,p,time,data);
        return WallDissipation(s,p,*current)-WallDissipation(s,p,*baseline);
    },face);
    // This is a local constrained native trace projection (endpoint/edge
    // preserving segment transforms, commuting tensor transforms in 3D).
    // Analytic chi is never projected. Only its physical dissipation change
    // is projected, with designated downstream face modes constrained to zero.
    for(int id:m_baselineTraceIds)trace[1].at(id)=0;
    return trace;
}

RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::ExternalLift(const Vector &state,
    double time,const RANSWallEvolutionData &data) const
{
    auto s=state;for(auto [row,value]:ExternalConstraints(time,data))
    {
        if(m_transitionEnabled && row>=offset[Epsilon()] && row<offset[InnerPressure()] && m_traceIds.count(row-offset[Epsilon()]))continue;
        s[row]=value;
    }
    return s;
}

RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::DirectionalAction(const Vector &state,
    const Vector &direction,const Vector &history,double alpha,double time,const RANSWallEvolutionData &data,double h) const
{
    Require(direction.size()==state.size() && h>0 && std::isfinite(h),"invalid directional action");
    auto plus=state,minus=state;
    for(int i=0;i<Size();++i){plus[i]+=h*direction[i];minus[i]-=h*direction[i];}
    const auto a=Residual(plus,history,alpha,time,data),b=Residual(minus,history,alpha,time,data);
    Vector result(Size());for(int i=0;i<Size();++i)result[i]=(a[i]-b[i])/(2*h);
    return result;
}
} // namespace Nektar
