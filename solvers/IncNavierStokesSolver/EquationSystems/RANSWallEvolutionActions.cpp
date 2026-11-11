///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <algorithm>
#include <cmath>
#include <map>

namespace Nektar
{
using namespace SolverUtils;
namespace
{
using Vector=RANSWallAlgebra::Vector;
using RANSWallAlgebra::Require;
void PressureDirection(const Vector &direction,int size,int differential)
{
    Require(static_cast<int>(direction.size())==size,"pressure direction has wrong size");
    for(int i=0;i<size;++i)
        Require(std::isfinite(direction[i]) && (i>=differential || direction[i]==0),
            "pressure action requires a finite pressure/gauge-only direction");
}
void Finite(const Vector &values)
{
    for(double value:values)Require(std::isfinite(value),"nonfinite pressure action");
}
}

RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::AffinePressureAction(
    const Vector &direction,const std::map<int,double> &constraints,int region) const
{
    PressureDirection(direction,Size(),DifferentialSize());
    Require(region>=-1 && region<=1,"invalid pressure action region");
    for(const auto &element:innerVelocity->Elements())
        Require(m_parameters.normalOrder>=element.degree,
            "exact state-free affine pressure action requires normal quadrature order at least the velocity degree");
    // These fixed physical rules include every native radial interface and h.
    // The affine stress term is p times a mapped cofactor/basis derivative;
    // the roughness term is a fixed tangential gradient times phi*J. No mass,
    // advection, constitutive stress or scalar source enters this assembly.
    std::vector<long double> owned(Size());
    std::map<RANSWallCoordinate,RANSVector> forceCache;
    auto integrate=[&](const VolumePoint &q,bool inner)
    {
        const int pressure=inner?InnerPressure():OuterPressure();
        const double value=q.p.Value(direction,offset[pressure]);
        RANSVector force{};
        if(inner && Distance(q.v.coordinate)<=m_parameters.geometry.roughness)
        {
            auto H=q.v.coordinate;H[2]=m_parameters.geometry.height;
            auto found=forceCache.find(H);
            if(found==forceCache.end())
            {
                const auto metric=m_mapping.At(H);
                const auto gradient=Point(OuterPressure(),H,1).Gradient(direction,offset[OuterPressure()]);
                double normal=0;for(int a=0;a<3;++a)normal+=metric.normal[a]*gradient[a];
                for(int a=0;a<3;++a)force[a]=gradient[a]-metric.normal[a]*normal;
                found=forceCache.emplace(H,force).first;
            }
            force=found->second;
        }
        for(int a=0;a<m_dimension;++a)for(int j=0;j<static_cast<int>(q.v.ids.size());++j)
        {
            const int aa=PhysicalAxis(a);
            owned[offset[a]+q.v.ids[j]]-=static_cast<long double>(q.weight)*
                (static_cast<long double>(value)*q.v.derivative[aa][j]+
                 static_cast<long double>(force[aa])*q.v.phi[j]);
        }
    };
    if(region!=1)for(const auto &q:m_bodyInner)integrate(q,true);
    if(region!=0)for(const auto &q:m_outer)integrate(q,false);
    Vector result(Size());for(int i=0;i<Size();++i)result[i]=static_cast<double>(owned[i]);
    m_comm->AllReduce(result,LibUtilities::ReduceSum);
    // The mean/gauge block is replicated. Add it only after owned assembly.
    if(!DevelopingFlat())for(int side=0;side<2;++side)if(region<0 || side==region)
        for(int i=0;i<static_cast<int>(m_pressureMean[side].size());++i)
        {
            const int row=offset[InnerPressure()+side]+i;
            result[row]+=direction[offset[Gauge()]]*m_pressureMean[side][i];
            result[offset[Gauge()]]+=direction[row]*m_pressureMean[side][i];
        }
    for(auto [row,value]:constraints){(void)value;result.at(row)=0.;}
    Finite(result);return result;
}

RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::PressureDirectionalAction(
    const Vector &state,const Vector &direction,double time,const RANSWallEvolutionData &data,
    bool imposeConstraints,int region) const
{
    PressureDirection(direction,Size(),DifferentialSize());
    Require(static_cast<int>(state.size())==Size(),"pressure action state has wrong size");
    const auto constraints=Constraints(state,time,data);
    if(AcceptedAWF())return AffinePressureAction(direction,imposeConstraints?constraints:std::map<int,double>{},region);
    auto s=m_transitionEnabled?ExternalLift(state,time,data):state;
    if(!m_transitionEnabled)for(auto [row,value]:constraints)s[row]=value;
    auto result=AffinePressureAction(direction,imposeConstraints?constraints:std::map<int,double>{},region);
    if(!m_transitionEnabled || m_parameters.chi==0)return result;

    auto chiDirection=[&](const RANSWallCoordinate &point)
    {
        auto H=point;H[2]=m_parameters.geometry.height;
        const auto input=TransitionInput(s,H,time,data);AWFTransitionInput derivative;
        const auto metric=m_mapping.At(H);
        const auto gradient=Point(OuterPressure(),H,1).Gradient(direction,offset[OuterPressure()]);
        double normal=0;for(int a=0;a<3;++a)normal+=metric.normal[a]*gradient[a];
        const double K0=SampleK(s,H),viscous=m_parameters.geometry.molecularNu/metric.height;
        const double scale2=K0/m_parameters.model.Q()+viscous*viscous;
        for(int a=0;a<3;++a)derivative.sensors.pressure[a]=metric.height*
            (gradient[a]-metric.normal[a]*normal)/scale2;
        return EvaluateEquilibriumAWFTransition(input,m_transition,derivative).derivative.chi;
    };
    struct RayDerivative
    {
        std::shared_ptr<EquilibriumAWFProfiles> profile;
        RANSCappedDestructionValues derivative;
        double chi=0,lambda=1;
    };
    std::map<RANSWallCoordinate,RayDerivative> cache;
    auto rayDerivative=[&](const RANSWallCoordinate &point)->const RayDerivative&
    {
        auto ground=point;ground[2]=0;
        const auto found=cache.find(ground);if(found!=cache.end())return found->second;
        RayDerivative value;value.chi=chiDirection(point);
        if(value.chi==0)return cache.emplace(ground,std::move(value)).first->second;
        value.profile=StateProfile(s,point,time,data);
        Require(value.profile->Reference().weight>0 && value.profile->Reference().cutoff>0,
            "nonzero pressure transition derivative requires an active qualified cutoff");
        const auto &profile=*value.profile;const double scale=m_mapping.At(point).scale;
        RANSCappedDestructionInput input;RANSCappedDestructionDirection derivative;
        input.sampledK=profile.K();
        for(const auto &q:profile.Quadrature())
        {
            auto chart=point;chart[2]=q.r/scale;
            const auto pre=profile.AtProduction(q.r,0);
            const double original=profile.Original().OriginalDissipation(q.r);
            const double corrected=profile.Reference().numerator/std::max(q.r,profile.Reference().cutoff);
            input.k.push_back(Point(K(),chart).Value(s,offset[K()]));
            input.preE.push_back(pre.dissipation);input.preD.push_back(pre.dissipation-pre.redistribution);
            input.metricWeights.push_back(q.weight*q.jacobian);
            derivative.preE.push_back(value.chi*(corrected-original));
            derivative.preD.push_back(value.chi*(pre.referenceProduction-original));
        }
        if(m_cappedDestruction)
        {
            value.lambda=EvaluateRANSCappedDestruction(input,m_destruction).lambda;
            value.derivative=RANSCappedDestructionDirectional(input,derivative,m_destruction);
        }
        else
        {
            value.derivative.destruction=std::move(derivative.preD);
            value.derivative.lambda=0;
        }
        return cache.emplace(ground,std::move(value)).first->second;
    };

    std::vector<long double> owned(Size());
    if(region!=1)for(const auto &ray:m_rays)
    {
        const auto &derivative=rayDerivative(ray.coordinate);
        if(derivative.chi==0)continue;
        const auto &rule=derivative.profile->Quadrature();
        Require(rule.size()==derivative.derivative.destruction.size(),"pressure source direction shape differs");
        const double scale=m_mapping.At(ray.coordinate).scale;
        for(std::size_t i=0;i<rule.size();++i)
        {
            auto point=ray.coordinate;point[2]=rule[i].r/scale;const auto q=Point(K(),point);
            if(Owner(q.geometry)!=m_comm->GetRank())continue;
            const long double weight=static_cast<long double>(ray.weight)*rule[i].weight*rule[i].jacobian;
            for(int j=0;j<static_cast<int>(q.ids.size());++j)
                owned[offset[K()]+q.ids[j]]+=weight*q.phi[j]*derivative.derivative.destruction[i];
        }
    }
    Vector nonlinear(Size());for(int i=0;i<Size();++i)nonlinear[i]=static_cast<double>(owned[i]);
    m_comm->AllReduce(nonlinear,LibUtilities::ReduceSum);
    for(int i=0;i<Size();++i)result[i]+=nonlinear[i];
    if(imposeConstraints)
    {
        for(auto [row,value]:constraints){(void)value;result[row]=0;}
        const auto face=[&](int axis,int side,const RANSWallGlobalField::Element &e)
            {return axis==2 && !side && std::abs(e.lower[2]-m_parameters.geometry.height)<1e-11;};
        const auto trace=fields[Epsilon()]->Boundary([&](const RANSWallCoordinate &point)
        {
            const auto &derivative=rayDerivative(point);if(derivative.chi==0)return 0.;
            const auto &profile=*derivative.profile;const double r=Distance(point);
            const double original=profile.Original().OriginalDissipation(r);
            const double corrected=profile.Reference().numerator/std::max(r,profile.Reference().cutoff);
            return derivative.derivative.lambda*profile.Dissipation(r)+
                derivative.lambda*derivative.chi*(corrected-original);
        },face);
        for(auto [id,value]:trace)
            result[offset[Epsilon()]+id]=m_baselineTraceIds.count(id)?0.:-value;
    }
    Finite(result);return result;
}
} // namespace Nektar
