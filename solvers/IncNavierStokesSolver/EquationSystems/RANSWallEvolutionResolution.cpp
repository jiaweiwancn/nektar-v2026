///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionResolution.h>
#include <LibUtilities/LinearAlgebra/Lapack.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>

namespace Nektar
{
using namespace SolverUtils;
using namespace RANSWallAlgebra;
namespace
{
using Wide=std::vector<long double>;
using Flux=std::array<double,3>;
struct Accumulator
{
    Vector sum,correction;Wide extended,absolute;Vector counts;
    explicit Accumulator(int n):sum(n),correction(n),extended(n),absolute(n),counts(n){}
    void Add(int row,double value,long double wide)
    {
        const double next=sum[row]+value;
        correction[row]+=std::abs(sum[row])>=std::abs(value)?
            (sum[row]-next)+value:(value-next)+sum[row];
        sum[row]=next;extended[row]+=wide;absolute[row]+=std::abs((long double)value);++counts[row];
    }
    Vector Compensated() const
    {
        Vector result(sum.size());for(int i=0;i<result.size();++i)result[i]=sum[i]+correction[i];return result;
    }
};
using Buckets=std::map<std::string,Accumulator>;
Accumulator &Bucket(Buckets &values,const std::string &name,int n)
{return values.try_emplace(name,n).first->second;}
void Add(Buckets &values,const std::string &name,int n,int offset,
    const RANSWallGlobalPoint &q,double value,const Flux &flux,double weight)
{
    auto &a=Bucket(values,name,n);
    for(int j=0;j<q.ids.size();++j)
    {
        // Match the baseline Add exactly: its initial value product is double,
        // flux products are long double, and each weighted increment is double.
        long double native=q.phi[j]*value;
        long double wide=(long double)q.phi[j]*value;
        for(int axis=0;axis<3;++axis)
        {native+=(long double)q.derivative[axis][j]*flux[axis];wide+=(long double)q.derivative[axis][j]*flux[axis];}
        a.Add(offset+q.ids[j],weight*static_cast<double>(native),(long double)weight*wide);
    }
}
Vector Cast(const Wide &values)
{Vector result;for(long double x:values){Require(std::isfinite(x),"nonfinite extended residual audit");result.push_back(x);}return result;}
Vector Difference(const Vector &a,const Vector &b)
{
    Require(a.size()==b.size(),"resolution vector shapes differ");Vector r(a.size());
    for(int i=0;i<r.size();++i){r[i]=a[i]-b[i];}
    return r;
}
class MassFactor
{
public:
    explicit MassFactor(const RANSWallAlgebra::DenseMatrix &a):n(a.size()),lu(n*n),pivots(n)
    {
        Require(n>0,"empty scalar free mass");
        for(int i=0;i<n;++i)for(int j=0;j<n;++j)lu[i+n*j]=a[i][j];
        int info=0;Lapack::Dgetrf(n,n,lu.data(),n,pivots.data(),info);Require(info==0,"singular scalar free mass");
    }
    Vector Solve(Vector rhs) const
    {auto pivot=pivots;int info=0;Lapack::Dgetrs('N',n,1,lu.data(),n,pivot.data(),rhs.data(),n,info);Require(info==0,"scalar mass-rate solve failed");return rhs;}
private:int n;Vector lu;std::vector<int> pivots;
};
}
const std::vector<std::string> &RANSWallResolutionSpatialTerms()
{
    static const std::vector<std::string> names{
        "inner_momentum_stress","outer_momentum_stress","conservative_advection",
        "k_diffusion","epsilon_diffusion","inner_k_production","inner_k_destruction",
        "outer_k_production","outer_k_destruction","epsilon_production",
        "epsilon_destruction","epsilon_yap","prescribed_scalar_source",
        "roughness_force","prescribed_body","continuity","gauge"};
    return names;
}

RANSWallEvolutionPatch::ResidualResolutionResult RANSWallEvolutionPatch::ResidualResolutionComponents(
    const Vector &state,const Vector &history,double alpha,double time,const RANSWallEvolutionData &data,bool anchored,
    bool originalOuterCounterfactual,bool analyticalReference) const
{
    Require(m_comm->GetSize()==1 && m_dimension==2 && !Mapped() && ((m_transitionEnabled && m_cappedDestruction) || AcceptedAWF()),
        "R3 residual-resolution audit requires serial flat2D Upstream+Capped");
    Require(state.size()==Size() && history.size()==Size() && alpha>=0 && m_eta==0,
        "R3 residual-resolution audit requires frozen unanchored valid state/history");
    (void)anchored;
    Require(!analyticalReference || (AcceptedAWF() && alpha==0),"analytical weak loads require the flat accepted neutral alpha=0 diagnostic");
    Require(!analyticalReference || (std::abs(m_referenceFlow[0]-1)<1e-14 && std::abs(m_referenceNormal[2]-1)<1e-14),
        "analytical weak loads currently support the declared +x flat neutral reference only");
    const double us=m_robustUstar,H=m_parameters.geometry.height,K0=m_parameters.model.Q()*us*us;
    const auto exact=analyticalReference?Profile(K0,{0,0,0},1):nullptr;
    double pressureJump=0.;
    if(analyticalReference)
    {
        RANSTensor G{};G[0][2]=us/(m_parameters.model.Kappa()*H);
        const auto c=EvaluateRANSClosure(m_parameters.model.Model(),K0,exact->Dissipation(H),G,H,m_parameters.model.Parameters(),0);
        pressureJump=c.reynoldsStress[2][2]-2*K0/3;
    }
    ResidualResolutionResult result;result.constraints=Constraints(state,time,data);
    result.effectiveState=ExternalLift(state,time,data);const auto &s=result.effectiveState;
    Refresh(s,time,&data);Buckets buckets;
    for(const auto &name:RANSWallResolutionSpatialTerms())Bucket(buckets,name,Size());
    auto load=[&](const std::string &name,int block,const RANSWallGlobalPoint &q,double value,const Flux &flux,double weight)
        {Add(buckets,name,Size(),offset[block],q,value,flux,weight);};
    auto evaluate=[&](const VolumePoint &q,bool inner,const EquilibriumAWFProfiles *profile,
        const RANSVector &force,const RANSCappedDestructionValues *normalization)
    {
        RANSTensor G{},stress{},cancelled{};RANSVector u{};
        for(int a=0;a<m_dimension;++a){u[PhysicalAxis(a)]=q.v.Value(s,offset[a]);G[PhysicalAxis(a)]=q.v.Gradient(s,offset[a]);}
        double k=q.v.Value(s,offset[K()]),pi=q.p.Value(s,offset[inner?InnerPressure():OuterPressure()]);
        auto gk=q.v.Gradient(s,offset[K()]);double dk=0,de=0,epsilon=0,sk=0,se=0;
        if(analyticalReference)
        {
            const double z=Distance(q.v.coordinate),UH=exact->Reference().matchingVelocity;
            u={};G={};gk={};
            u[0]=inner?exact->Original().Velocity(z,{UH,0},{0,0})[0]:UH+us/m_parameters.model.Kappa()*std::log(z/H);
            G[0][2]=inner?exact->Original().VelocityGradient(z,{UH,0},{0,0})[0]:us/(m_parameters.model.Kappa()*z);
            k=inner?exact->ReferenceK(z):K0;gk[2]=inner?exact->ReferenceKGradient(z):0.;
            pi=inner?pressureJump*(upper[2]-H)/upper[2]:-pressureJump*H/upper[2];
            if(inner)profile=exact.get();
        }
        if(inner)
        {
            const double nt=profile->Original().EddyViscosity(Distance(q.v.coordinate));
            const auto law=InnerStress(k,nt,G);stress=RANSWallInnerMeanStress(law,m_parameters.geometry.molecularNu,pi);
            dk=m_parameters.geometry.molecularNu+nt/m_parameters.model.Parameters().sigmaK;
            const auto actual=WallSource(s,q.v.coordinate,law.production,*profile,normalization);
            double destruction=actual.dissipation-actual.redistribution;
            if(m_cappedDestruction)
            {
                const auto pre=profile->AtProduction(Distance(q.v.coordinate),law.production);
                destruction=EvaluateRANSCappedDestructionPoint(k,profile->K(),pre.dissipation,
                    pre.dissipation-pre.redistribution,normalization->b,normalization->m,m_destruction).destruction;
            }
            sk=actual.kSource+data.innerKSource(q.v.coordinate,time);
            load("inner_k_production",K(),q.v,-law.production,{},q.weight);
            load("inner_k_destruction",K(),q.v,destruction,{},q.weight);
            // These auxiliary E and -C loads are NEVER in the physical sum.
            load("diagnostic_inner_Eplus",K(),q.v,actual.dissipation,{},q.weight);
            load("diagnostic_inner_minus_Cplus",K(),q.v,-actual.redistribution,{},q.weight);
            load("prescribed_scalar_source",K(),q.v,-data.innerKSource(q.v.coordinate,time),{},q.weight);
        }
        else
        {
            epsilon=analyticalReference?us*us*us/(m_parameters.model.Kappa()*Distance(q.v.coordinate)):q.e.Value(s,offset[Epsilon()]);const auto c=originalOuterCounterfactual?
                EvaluateRANSClosure(m_parameters.model.Model(),k,epsilon,G,Distance(q.v.coordinate),m_parameters.model.Parameters(),0):
                OuterClosure(k,epsilon,G,q.v.coordinate);
            dk=c.kDiffusion;de=c.epsilonDiffusion;sk=c.kSource+data.outerKSource(q.v.coordinate,time);
            se=c.epsilonSource+data.outerEpsilonSource(q.v.coordinate,time);
            const double distance=Distance(q.v.coordinate),kap=m_parameters.model.Kappa();
            const double gamma=originalOuterCounterfactual?0:m_regularizationGamma+m_regularizationScale*kap*m_robustUstar*distance;
            const double shear=m_robustUstar/(kap*distance);long double stableProduction=0;
            for(int a=0;a<3;++a)for(int b=0;b<3;++b)
            {
                stress[a][b]=-c.reynoldsStress[a][b]+(a==b?2*k/3-pi:0);
                const double sr=shear*(m_referenceFlow[a]*m_referenceNormal[b]+m_referenceNormal[a]*m_referenceFlow[b]);
                const long double B=(long double)c.eddyViscosity*(G[a][b]+G[b][a])-c.quadraticStress[a][b]-c.cubicStress[a][b]+
                    (long double)gamma*(G[a][b]+G[b][a]-sr);
                cancelled[a][b]=B-(a==b?pi:0);
                stableProduction+=(B-(a==b?(long double)2*k/3:0))*G[a][b];
            }
            load("outer_k_production",K(),q.v,-c.production,{},q.weight);
            load("outer_k_destruction",K(),q.v,epsilon,{},q.weight);
            load("epsilon_production",Epsilon(),q.e,-c.epsilonProduction,{},q.weight);
            load("epsilon_destruction",Epsilon(),q.e,c.epsilonDestruction,{},q.weight);
            load("epsilon_yap",Epsilon(),q.e,-c.yap,{},q.weight);
            load("prescribed_scalar_source",K(),q.v,-data.outerKSource(q.v.coordinate,time),{},q.weight);
            load("prescribed_scalar_source",Epsilon(),q.e,-data.outerEpsilonSource(q.v.coordinate,time),{},q.weight);
            load("diagnostic_outer_cancelled_k_production",K(),q.v,-static_cast<double>(stableProduction),{},q.weight);
            load("diagnostic_outer_cancelled_epsilon_production",Epsilon(),q.e,
                -static_cast<double>((long double)m_parameters.model.Parameters().cEpsilon1*epsilon/k*stableProduction),{},q.weight);
        }
        for(int a=0;a<m_dimension;++a)
        {
            const int aa=PhysicalAxis(a);Flux advection{},nativeFlux{};
            for(int b=0;b<3;++b){advection[b]=-u[aa]*u[b];nativeFlux[b]=stress[aa][b]-u[aa]*u[b];}
            const double rough=inner && Distance(q.v.coordinate)<=m_parameters.geometry.roughness?force[aa]:0;
            load(inner?"inner_momentum_stress":"outer_momentum_stress",a,q.v,0,stress[aa],q.weight);
            if(!inner)load("diagnostic_outer_cancelled_momentum",a,q.v,0,cancelled[aa],q.weight);
            load("conservative_advection",a,q.v,0,advection,q.weight);
            load("roughness_force",a,q.v,-rough,{},q.weight);
            load("legacy_spatial",a,q.v,-rough,nativeFlux,q.weight);
        }
        Flux diffusion{},advection{},nativeK{};
        for(int b=0;b<3;++b){diffusion[b]=dk*gk[b];advection[b]=-u[b]*k;nativeK[b]=dk*gk[b]-u[b]*k;}
        load("k_diffusion",K(),q.v,0,diffusion,q.weight);
        load("conservative_advection",K(),q.v,0,advection,q.weight);
        load("legacy_spatial",K(),q.v,-sk,nativeK,q.weight);
        const double divergence=G[0][0]+G[1][1]+G[2][2];
        load("continuity",inner?InnerPressure():OuterPressure(),q.p,divergence,{},q.weight);
        load("legacy_spatial",inner?InnerPressure():OuterPressure(),q.p,divergence,{},q.weight);
        if(!inner)
        {
            auto ge=q.e.Gradient(s,offset[Epsilon()]);Flux ed{},ea{},ef{};
            if(analyticalReference){ge={};ge[2]=-epsilon/Distance(q.v.coordinate);}
            for(int b=0;b<3;++b){ed[b]=de*ge[b];ea[b]=-u[b]*epsilon;ef[b]=de*ge[b]-u[b]*epsilon;}
            load("epsilon_diffusion",Epsilon(),q.e,0,ed,q.weight);
            load("conservative_advection",Epsilon(),q.e,0,ea,q.weight);
            load("legacy_spatial",Epsilon(),q.e,-se,ef,q.weight);
        }
    };
    for(const auto &ray:m_rays)
    {const auto force=analyticalReference?RANSVector{}:CartesianForcing(s,ray.coordinate,time,data,false);for(const auto &q:ray.points)evaluate(q,true,ray.profile.get(),force,&ray.destruction);}
    for(const auto &q:m_outer)evaluate(q,false,nullptr,{},nullptr);
    const auto prescribed=OwnedBodyLoad(time,data);
    for(int row=0;row<Size();++row)
    {Bucket(buckets,"prescribed_body",Size()).Add(row,-prescribed[row],-prescribed[row]);Bucket(buckets,"legacy_spatial",Size()).Add(row,-prescribed[row],-prescribed[row]);}
    for(int side=0;side<2;++side)for(int j=0;j<m_pressureMean[side].size();++j)
    {
        const int row=offset[InnerPressure()+side]+j;const double w=m_pressureMean[side][j];
        for(const std::string name:{"gauge","legacy_spatial"})
        {Bucket(buckets,name,Size()).Add(row,s[offset[Gauge()]]*w,(long double)s[offset[Gauge()]]*w);
         Bucket(buckets,name,Size()).Add(offset[Gauge()],s[row]*w,(long double)s[row]*w);}
    }
    auto &mass=Bucket(buckets,"time_mass",Size());
    auto &massWideDelta=Bucket(buckets,"diagnostic_time_mass_extended_subtraction",Size());
    for(int block=0;block<InnerPressure();++block)
    {
        Vector delta(size[block]);for(int j=0;j<size[block];++j)delta[j]=s[offset[block]+j]-history[offset[block]+j];
        const auto baseline=Multiply(m_mass[block],delta);
        for(int i=0;i<size[block];++i)
        {
            long double wide=0,wideDelta=0;
            for(int j=0;j<size[block];++j)
            {wide+=(long double)m_mass[block][i][j]*delta[j];
             wideDelta+=(long double)m_mass[block][i][j]*((long double)s[offset[block]+j]-history[offset[block]+j]);}
            mass.Add(offset[block]+i,alpha*baseline[i],(long double)alpha*wide);
            massWideDelta.Add(offset[block]+i,alpha*baseline[i],(long double)alpha*wideDelta);
        }
    }
    for(auto &[name,values]:buckets)
    {
        result.terms[name]=values.sum;result.compensated[name]=values.Compensated();result.extended[name]=values.extended;
        result.absoluteAccumulation[name]=Vector(Size());result.additionCounts[name]=values.counts;
        for(int i=0;i<Size();++i)
        {
            // Enclose the positive long-double sum of frozen binary64
            // increments; this does not bound errors in their physical inputs.
            const long double u=std::numeric_limits<long double>::epsilon()/2;
            const long double denominator=1-values.counts[i]*u;
            Require(denominator>0,"invalid resolution accumulation bound");
            const double upper=values.absolute[i]==0?0:static_cast<double>(values.absolute[i]/denominator);
            result.absoluteAccumulation[name][i]=upper==0?0:std::nextafter(std::nextafter(upper,
                std::numeric_limits<double>::infinity()),std::numeric_limits<double>::infinity());
        }
    }
    auto total=[&](const std::string &key,const std::vector<std::string> &names)
    {
        result.terms[key]=Vector(Size());result.compensated[key]=Vector(Size());result.extended[key]=Wide(Size());
        for(const auto &name:names)for(int i=0;i<Size();++i)
        {result.terms[key][i]+=result.terms.at(name)[i];result.compensated[key][i]+=result.compensated.at(name)[i];result.extended[key][i]+=result.extended.at(name)[i];}
    };
    total("spatial_sum",RANSWallResolutionSpatialTerms());
    total("total_unmasked",{"spatial_sum","time_mass"});
    result.extended["total_unmasked_cancelled"]=result.extended["total_unmasked"];
    for(int i=0;i<Size();++i)for(const auto &pair:std::vector<std::pair<std::string,std::string>>{
        {"diagnostic_outer_cancelled_momentum","outer_momentum_stress"},
        {"diagnostic_outer_cancelled_k_production","outer_k_production"},
        {"diagnostic_outer_cancelled_epsilon_production","epsilon_production"}})
        result.extended["total_unmasked_cancelled"][i]+=result.extended.at(pair.first)[i]-result.extended.at(pair.second)[i];
    result.terms["legacy_stage"]=result.terms["legacy_spatial"];
    result.compensated["legacy_stage"]=result.compensated["legacy_spatial"];
    result.extended["legacy_stage"]=result.extended["legacy_spatial"];
    const auto compensatedMass=mass.Compensated();
    for(int i=0;i<Size();++i)
    {result.terms["legacy_stage"][i]+=mass.sum[i];result.compensated["legacy_stage"][i]+=compensatedMass[i];result.extended["legacy_stage"][i]+=mass.extended[i];}
    for(const std::string key:{"total_constrained","total_constrained_cancelled"})
    {
        const std::string parent=key=="total_constrained"?"total_unmasked":"total_unmasked_cancelled";
        result.extended[key]=result.extended[parent];
        if(key=="total_constrained"){result.terms[key]=result.terms[parent];result.compensated[key]=result.compensated[parent];}
        for(auto [row,value]:result.constraints)
        {
            result.extended[key][row]=(long double)state[row]-value;
            if(key=="total_constrained"){result.terms[key][row]=state[row]-value;result.compensated[key][row]=state[row]-value;}
        }
    }
    result.terms["constraint_replacement"]=Difference(result.terms["total_constrained"],result.terms["total_unmasked"]);
    result.extended["constraint_replacement"]=Wide(Size());
    for(int i=0;i<Size();++i)result.extended["constraint_replacement"][i]=result.extended["total_constrained"][i]-result.extended["total_unmasked"][i];
    result.scope="Diagnostic arithmetic at fixed state/history/actual quadrature; net destruction appears once. Extended products retain native rounded physical primitives and stored mass entries. No roundoff-floor, positivity, trajectory or root claim.";
    return result;
}

RANSWallEvolutionPatch::ResidualResolutionResult RANSWallEvolutionPatch::ResidualResolution(
    const Vector &state,const Vector &history,double alpha,double time,const RANSWallEvolutionData &data,bool anchored) const
{
    auto result=ResidualResolutionComponents(state,history,alpha,time,data,anchored);
    result.native["spatial_unmasked"]=Residual(state,history,0,time,data,anchored,false);
    result.native["stage_unmasked"]=Residual(state,history,alpha,time,data,anchored,false);
    result.native["stage_constrained"]=Residual(state,history,alpha,time,data,anchored,true);
    result.native["external_lift_delta"]=Difference(result.effectiveState,state);
    result.norms["native_stage"]=Norms(result.native["stage_constrained"],result.constraints,m_robustUstar);
    for(const std::string name:{"total_constrained","legacy_stage"})
    {
        auto values=result.terms[name];for(auto [row,value]:result.constraints)values[row]=state[row]-value;
        result.norms[name]=Norms(values,result.constraints,m_robustUstar);
        result.norms["native_minus_"+name]=Norms(Difference(result.native["stage_constrained"],values),result.constraints,m_robustUstar);
    }
    for(const std::string name:{"total_constrained","total_constrained_cancelled"})
    {
        const auto values=Cast(result.extended[name]);
        result.norms["extended_"+name]=Norms(values,result.constraints,m_robustUstar);
        result.norms["native_minus_extended_"+name]=Norms(Difference(result.native["stage_constrained"],values),result.constraints,m_robustUstar);
    }
    for(const auto &name:RANSWallResolutionSpatialTerms())result.norms["term_"+name]=Norms(result.terms[name],result.constraints,m_robustUstar);
    return result;
}

RANSWallScalarMassRate RANSWallResolutionMassRate(const RANSWallEvolutionPatch &patch,
    const RANSWallEvolutionPatch::ResidualResolutionResult &resolution,int block,
    const RANSWallCoordinate &contact,const Vector &boundaryRate)
{
    Require(block==patch.K() || block==patch.Epsilon(),"contact rate requires a scalar block");
    const int n=patch.size[block],offset=patch.offset[block];const auto mass=patch.fields[block]->Mass();
    RANSWallScalarMassRate result;result.boundaryRate=boundaryRate.empty()?Vector(n):boundaryRate;
    Require(result.boundaryRate.size()==n,"scalar boundary-rate shape differs");
    std::vector<int> free,fixed;
    for(int i=0;i<n;++i)(resolution.constraints.count(offset+i)?fixed:free).push_back(i);
    for(int i:free)Require(result.boundaryRate[i]==0,"boundary-rate input has a nonzero free entry");
    for(int i:free){result.freeRows.push_back(i);}
    for(int i:fixed){result.fixedRows.push_back(i);}
    RANSWallAlgebra::DenseMatrix ff(free.size(),Vector(free.size()));
    for(int i=0;i<free.size();++i)for(int j=0;j<free.size();++j)ff[i][j]=mass[free[i]][free[j]];
    const MassFactor factor(ff);const auto &spatial=resolution.native.at("spatial_unmasked");
    result.spatialLoad.assign(spatial.begin()+offset,spatial.begin()+offset+n);
    auto solve=[&](const Vector &load,bool boundary)
    {
        Vector rhs(free.size());for(int i=0;i<free.size();++i)
        {rhs[i]=-load[free[i]];if(boundary)for(int j:fixed)rhs[i]-=mass[free[i]][j]*result.boundaryRate[j];}
        const auto value=factor.Solve(rhs);Vector rate(n);for(int i=0;i<free.size();++i)rate[free[i]]=value[i];
        if(boundary)for(int i:fixed){rate[i]=result.boundaryRate[i];}
        return rate;
    };
    result.rate=solve(result.spatialLoad,true);
    for(const auto &name:RANSWallResolutionSpatialTerms())
    {const auto &full=resolution.terms.at(name);Vector load(full.begin()+offset,full.begin()+offset+n);result.componentRates[name]=solve(load,false);}
    result.componentRates["boundary_rate_mass"]=solve(Vector(n),true);
    const auto applied=Multiply(mass,result.rate);for(int i:free)result.freeMassResidual.push_back(applied[i]+result.spatialLoad[i]);
    const auto row=patch.Point(block,contact,block==patch.Epsilon() && std::abs(contact[2]-patch.Parameters().geometry.height)<1e-11?1:0);
    result.contactValue=row.Value(resolution.effectiveState,offset);result.contactRate=row.Value(result.rate);
    for(const auto &[name,rate]:result.componentRates)result.contactComponents[name]=row.Value(rate);
    return result;
}
} // namespace Nektar
