///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallTiming.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallManufactured.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionResolution.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallForwardDiagnostics.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <algorithm>
#include <cmath>
#include <limits>

namespace Nektar
{
using namespace SolverUtils;
using RANSWallAlgebra::Require;
namespace
{
double MinimumEigenvalue(RANSTensor a)
{
    return RANSWallSolveJson::Eigenvalues(a)[0];
}
}
RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::NeutralAccuracy(const Vector &s,double time,const RANSWallEvolutionData &data) const
{
    RANSWallTiming::Scope timing("audit.neutral_accuracy");
    Require(NewModels() && m_dimension==2 && !Mapped() && m_comm->GetSize()==1,"neutral accuracy requires serial flat 2D new models");
    Require(m_referenceFlow==RANSVector{1,0,0} && m_referenceNormal==RANSVector{0,0,1},
        "NeutralV1 accuracy observer requires canonical positive-x reference direction");
    const double us=m_robustUstar,H=m_parameters.geometry.height,Z=upper[2],K0=m_parameters.model.Q()*us*us;
    const double A=us*us*us/m_parameters.model.Kappa();const auto profile=Profile(K0);
    const double UH=profile->Reference().matchingVelocity;
    RANSTensor G{};G[0][2]=us/(m_parameters.model.Kappa()*H);
    const auto c=OuterClosure(K0,A/H,G,{0,0,H});
    const double jump=c.reynoldsStress[2][2]-2*K0/3;
    const std::array<double,2> pi{jump*(Z-H)/Z,-jump*H/Z};
    if(m_accuracyPoints.empty())for(const auto &q:Quadrature(0))
    {
        AccuracyPoint p;p.velocity=q;p.outer=q.coordinate[2]>H;
        p.pressure=Point(InnerPressure()+p.outer,q.coordinate);
        if(p.outer)p.epsilon=Point(Epsilon(),q.coordinate);
        const double z=q.coordinate[2];
        p.exact={p.outer?UH+us/m_parameters.model.Kappa()*std::log(z/H):profile->Original().Velocity(z,{UH,0},{0,0})[0],
            0.,p.outer || !AcceptedAWF()?K0:profile->ReferenceK(z),p.outer?A/z:0.,pi[p.outer]};
        p.initial={q.Value(m_reference,offset[0]),q.Value(m_reference,offset[1]),q.Value(m_reference,offset[K()]),
            p.outer?p.epsilon.Value(m_reference,offset[Epsilon()]):0.,p.pressure.Value(m_reference,offset[InnerPressure()+p.outer])};
        m_accuracyPoints.push_back(std::move(p));
    }
    const std::array<double,5> scale{us,us,K0,A/H,us*us};
    std::array<long double,5> e{},d{},volume{};long double pe=0,pd=0;
    for(const auto &p:m_accuracyPoints)
    {
        const double v=p.pressure.Value(s,offset[InnerPressure()+p.outer]);
        pe+=p.velocity.weight*(v-p.exact[4]);pd+=p.velocity.weight*(v-p.initial[4]);
        volume[4]+=p.velocity.weight;
    }
    const double eg=pe/volume[4],dg=pd/volume[4];volume[4]=0;
    for(const auto &p:m_accuracyPoints)
    {
        const auto &q=p.velocity;
        const std::array<double,5> v{q.Value(s,offset[0]),q.Value(s,offset[1]),q.Value(s,offset[K()]),
            p.outer?p.epsilon.Value(s,offset[Epsilon()]):0.,p.pressure.Value(s,offset[InnerPressure()+p.outer])};
        for(int j=0;j<5;++j)if(j!=3 || p.outer)
        {
            const double a=(v[j]-p.exact[j]-(j==4?eg:0.))/scale[j],b=(v[j]-p.initial[j]-(j==4?dg:0.))/scale[j];
            e[j]+=q.weight*a*a;d[j]+=q.weight*b*b;volume[j]+=q.weight;
        }
    }
    Vector result(12);for(int j=0;j<5;++j){result[j]=std::sqrt(e[j]/volume[j]);result[5+j]=std::sqrt(d[j]/volume[j]);}
    RANSVector wall{};double area=0;
    for(const auto &f:Surface(0,2,0.,-1,8))
    {
        const auto &q=f.test;RANSTensor g{};for(int a=0;a<2;++a)g[PhysicalAxis(a)]=q.Gradient(s,offset[a]);
        const auto live=StateProfile(s,q.coordinate,time,data);
        const auto t=RANSWallInnerMeanStress(InnerStress(q.Value(s,offset[K()]),live->Original().EddyViscosity(0.),g),m_parameters.geometry.molecularNu,
            Point(InnerPressure(),q.coordinate).Value(s,offset[InnerPressure()])-eg);
        const double w=f.weight*std::hypot(f.area[0],f.area[1],f.area[2]);area+=w;
        for(int a=0;a<3;++a)wall[a]+=w*(t[a][2]-(a==0?us*us:a==2?-pi[0]:0.));
    }
    for(double v:wall)result[10]+=v*v;
    result[10]=std::sqrt(result[10])/area/(us*us);result[11]=eg;
    Require(std::all_of(result.begin(),result.end(),[](double v){return std::isfinite(v);}),"nonfinite neutral accuracy measurement");
    return result;
}
RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::ManufacturedAccuracy(const Vector &s,const Vector &initial,
    double time,const RANSWallEvolutionData &data) const
{
    RANSWallTiming::Scope timing("audit.manufactured_accuracy");
    Require(CurvedManufactured() && initial.size()==s.size(),"CT6 accuracy requires its manufactured route and actual initial vector");
    const RANSWallManufactured exact(*this);const double K0=m_parameters.model.Q()*.16;
    const std::array<double,5> scale{.2,.2,K0,.4*.4*.4/(m_parameters.model.Kappa()*m_parameters.geometry.height),.1};
    std::array<long double,5> error{},drift{},volume{};
    for(const auto &v:Quadrature(0,6))
    {
        const bool outer=v.coordinate[2]>m_parameters.geometry.height;const auto truth=exact.At(v.coordinate,time,outer);
        const auto ep=outer?GeometryPoint(Epsilon(),v.geometry,v.coordinate):RANSWallGlobalPoint{};
        const auto pressure=GeometryPoint(InnerPressure()+outer,v.geometry,v.coordinate);
        const std::array<double,5> actual{v.Value(s,offset[0]),v.Value(s,offset[1]),v.Value(s,offset[K()]),
            outer?ep.Value(s,offset[Epsilon()]):0.,pressure.Value(s,offset[InnerPressure()+outer])};
        const std::array<double,5> start{v.Value(initial,offset[0]),v.Value(initial,offset[1]),v.Value(initial,offset[K()]),
            outer?ep.Value(initial,offset[Epsilon()]):0.,pressure.Value(initial,offset[InnerPressure()+outer])};
        for(int b=0;b<5;++b)if(b!=3 || outer)
        {
            const double e=(actual[b]-truth.value[b])/scale[b],d=(actual[b]-start[b])/scale[b];
            error[b]+=v.weight*e*e;drift[b]+=v.weight*d*d;volume[b]+=v.weight;
        }
    }
    Vector result(15);for(int b=0;b<5;++b){result[b]=std::sqrt(error[b]/volume[b]);result[b+5]=std::sqrt(drift[b]/volume[b]);}
    long double ground=0,groundArea=0,jump=0,innerError=0,outerError=0,interfaceArea=0;
    const auto stress=[&](const RANSWallCoordinate &q,int side)
    {
        const int traceSide=q[2]==m_parameters.geometry.height?side:0;
        const auto v=Point(0,q,traceSide);RANSTensor G{},T{};for(int a=0;a<2;++a)G[PhysicalAxis(a)]=v.Gradient(s,offset[a]);
        const double k=v.Value(s,offset[K()]),pi=Point(side<0?InnerPressure():OuterPressure(),q,traceSide).Value(s,offset[side<0?InnerPressure():OuterPressure()]);
        if(side<0)
        {const auto f=StateProfile(s,q,time,data);T=RANSWallInnerMeanStress(InnerStress(k,f->Original().EddyViscosity(Distance(q)),G),m_parameters.geometry.molecularNu,pi);}
        else
        {
            const double e=Point(Epsilon(),q,1).Value(s,offset[Epsilon()]);const auto law=OuterClosure(k,e,G,q);
            for(int i=0;i<3;++i)for(int j=0;j<3;++j)T[i][j]=-law.reynoldsStress[i][j]+(i==j?2*k/3-pi:0.);
        }
        return T;
    };
    for(const auto &face:Surface(0,2,0.,-1,12))
    {
        const auto &q=face.test.coordinate;const auto actual=stress(q,-1),truth=exact.At(q,time,false).stress;
        const double measure=std::hypot(face.area[0],face.area[1],face.area[2]),w=face.weight*measure;groundArea+=w;
        for(int i=0;i<3;++i){double d=0;for(int j=0;j<3;++j)d+=(actual[i][j]-truth[i][j])*face.area[j]/measure;ground+=w*d*d;}
    }
    for(const auto &face:Surface(0,2,m_parameters.geometry.height,1,12))
    {
        const auto &q=face.test.coordinate;const auto inner=stress(q,-1),outer=stress(q,1),truth=exact.At(q,time,false).stress;
        const double measure=std::hypot(face.area[0],face.area[1],face.area[2]),w=face.weight*measure;interfaceArea+=w;
        for(int i=0;i<3;++i)
        {
            double d=0,di=0,dout=0;
            for(int j=0;j<3;++j){const double n=face.area[j]/measure;d+=(inner[i][j]-outer[i][j])*n;di+=(inner[i][j]-truth[i][j])*n;dout+=(outer[i][j]-truth[i][j])*n;}
            jump+=w*d*d;innerError+=w*di*di;outerError+=w*dout*dout;
        }
    }
    result[10]=std::sqrt(ground/groundArea)/.16;result[11]=std::sqrt(jump/interfaceArea)/.16;
    result[12]=std::sqrt(innerError/interfaceArea)/.16;result[13]=std::sqrt(outerError/interfaceArea)/.16;
    result[14]=6;return result;
}
RANSWallEvolutionPatch::Diagnostic RANSWallEvolutionPatch::Audit(const Vector &s,const Vector &history,
    double alpha,double time,const RANSWallEvolutionData &data,double ustar,bool constitutiveReconstruction) const
{
    RANSWallTiming::Scope timing("audit.physical");
    if(ConstitutiveKInventory() && !m_assertKCoverage)
    {
        const auto domain=ScalarDomain(s,time,data);
        Require(domain.admissible,"audit outside ConstitutiveV2 inventory: "+domain.failure);
    }
    CoverageScope coverage(m_assertKCoverage,ConstitutiveKInventory());
    Diagnostic result;const double sentinel=std::numeric_limits<double>::max();
    Vector range{0,sentinel,-sentinel,sentinel,-sentinel,sentinel,-sentinel,0};
    Vector zeroInnerCovariance{0,sentinel,-sentinel};
    const double K0=m_parameters.model.Q()*ustar*ustar;
    bool admissible=true;Refresh(s,time,&data);
    for(const std::string name:{"source_volume_inner","source_volume_outer","physical_sample_K","physical_sample_outer",
        "prescribed_source_inner","mass_volume","interface_inner","interface_outer","reconstruction_dense","reconstruction_native"})
        result[name]={0,sentinel,sentinel,sentinel,sentinel,0};
    // Every rank evaluates replicated diagnostic sets. Owned source sets are
    // reduced below; replicated sets are counted only on rank zero.
    auto sample=[&](const std::string &name,const RANSWallGlobalPoint &v,bool inner,bool constitutive,bool owned)
    {
        if(!owned && m_comm->GetRank()!=0)return;
        constitutive=constitutive || constitutiveReconstruction;
        auto &r=result[name];if(r.empty())r={0,sentinel,sentinel,sentinel,sentinel,0};
        if(m_assertKCoverage)Require(m_coveredKLocations.count(v.coordinate),"unprotected constitutive k evaluation: audit sample");
        ++r[0];const double k=v.Value(s,offset[K()]);r[1]=std::min(r[1],k);
        const double epsilon=inner?(constitutive?DependentEpsilon(s,v.coordinate,time,data):sentinel):
            GeometryPoint(Epsilon(),v.geometry,v.coordinate).Value(s,offset[Epsilon()]);
        r[2]=std::min(r[2],epsilon);
        if(!std::isfinite(k) || k<0 || (!inner && !(k>0)) || !std::isfinite(epsilon) || !(epsilon>0))++r[5];
        // The selected inner law is defined at exact local k=0.  Evaluate
        // its covariance there even when a collapsed source trace is itself
        // invalid for the outer high-Re equations; do not silently skip it.
        const bool stressState=std::isfinite(k) && (inner?k>=0:
            (k>0 && epsilon>0 && std::isfinite(epsilon)));
        if(constitutive && stressState)
        {
            RANSTensor G{},covariance{};for(int a=0;a<m_dimension;++a)G[PhysicalAxis(a)]=v.Gradient(s,offset[a]);
            if(inner)
            {
                const auto profile=StateProfile(s,v.coordinate,time,data);
                const double nt=profile->Original().EddyViscosity(Distance(v.coordinate));
                covariance=InnerStress(k,nt,G).reynoldsStress;
                if(constitutiveReconstruction)
                {range[5]=std::min(range[5],profile->K()/K0);range[6]=std::max(range[6],profile->K()/K0);}
            }
            else
            {
                const auto c=NewModels()?RANSWallDiagnostics::Observe(m_parameters.model,k,epsilon,G,Distance(v.coordinate)).closure:
                    OuterClosure(k,epsilon,G,v.coordinate);
                covariance=c.reynoldsStress;
                if(constitutiveReconstruction)
                {
                    const double g0=ustar/(m_parameters.model.Kappa()*Distance(v.coordinate));
                    const double e0=ustar*ustar*ustar/(m_parameters.model.Kappa()*Distance(v.coordinate));
                    for(int a=0;a<3;++a)for(int b=0;b<3;++b)
                        range[0]=std::max(range[0],std::abs(G[a][b]/g0-m_referenceFlow[a]*m_referenceNormal[b]));
                    range[1]=std::min(range[1],k/K0);range[2]=std::max(range[2],k/K0);
                    range[3]=std::min(range[3],epsilon/e0);range[4]=std::max(range[4],epsilon/e0);
                    range[7]=std::max(range[7],std::abs(G[0][0]+G[1][1]+G[2][2])/g0);
                }
                if(!NewModels())
                {
                    const double h=1e-5*std::max(1.,std::abs(G[0][2]));auto plus=G,minus=G;plus[1][0]+=h;minus[1][0]-=h;
                    const auto cp=OuterClosure(k,epsilon,plus,v.coordinate);
                    const auto cm=OuterClosure(k,epsilon,minus,v.coordinate);
                    r[4]=std::min(r[4],-(cp.reynoldsStress[1][0]-cm.reynoldsStress[1][0])/(2*h));
                }
            }
            const double minimum=MinimumEigenvalue(covariance);
            r[3]=std::min(r[3],minimum);
            if(inner && k==0)
            {
                ++zeroInnerCovariance[0];
                zeroInnerCovariance[1]=std::min(zeroInnerCovariance[1],minimum);
                zeroInnerCovariance[2]=std::max(zeroInnerCovariance[2],minimum);
            }
        }
    };
    for(const auto &ray:m_rays)
    {
        for(const auto &q:ray.points)sample("source_volume_inner",q.v,true,true,true);
        if(!ColumnAverageK())sample("physical_sample_K",ray.sample,true,true,false);
        else {auto &r=result["physical_sample_K"];++r[0];r[1]=std::min(r[1],SampleK(s,ray.coordinate));}
        if(m_transitionEnabled)
        {
            const double scale=m_mapping.At(ray.coordinate).scale;
            for(double distance:{m_sensorA,m_sensorB})
            {auto p=ray.coordinate;p[2]=distance/scale;sample("physical_sample_outer",Point(0,p),false,true,false);}
        }
    }
    for(const auto &q:m_outer)sample("source_volume_outer",q.v,false,true,true);
    for(const auto &q:m_bodyInner)sample("prescribed_source_inner",q.v,true,false,true);
    for(const auto &q:Quadrature(0))sample("mass_volume",q,q.coordinate[2]<m_parameters.geometry.height,false,false);
    for(const auto &f:Surface(0,2,m_parameters.geometry.height,1,8))
    {
        sample("interface_inner",Point(0,f.test.coordinate,-1),true,true,false);
        sample("interface_outer",Point(0,f.test.coordinate,1),false,true,false);
    }
    // Include actual native GLL output/reconstruction locations and a denser
    // equispaced set, independently of source Gauss points. This remains a
    // finite sample, even when every sampled value is positive.
    for(const auto &p:ReconstructionLocations(0,constitutiveReconstruction))
        sample(p.inventory,p.row,!p.strict,false,false);
    for(auto &[name,r]:result)
    {
        (void)name;for(int i:{0,5})m_comm->AllReduce(r[i],LibUtilities::ReduceSum);
        for(int i=1;i<=4;++i)m_comm->AllReduce(r[i],LibUtilities::ReduceMin);
        if(r[5])admissible=false;
        // Finite JSON placeholder for quantities not evaluated on this set.
        for(int i=1;i<=4;++i)if(r[i]==sentinel)r[i]=0;
    }
    result["admissible"]={admissible?1.:0.};
    // r[4] is a legacy tangent observation only. The new closures are forward
    // qualified: an unavailable probe must neither invalidate a physical A95
    // state nor be interpreted as a positive ellipticity certificate.
    result["outer_tangent_diagnostic_available"]={NewModels()?0.:1.};
    if(constitutiveReconstruction)
    {
        for(int i:{0,2,4,6,7})m_comm->AllReduce(range[i],LibUtilities::ReduceMax);
        for(int i:{1,3,5})m_comm->AllReduce(range[i],LibUtilities::ReduceMin);
        result["outer_constitutive_range"]=range;
        m_comm->AllReduce(zeroInnerCovariance[0],LibUtilities::ReduceSum);
        m_comm->AllReduce(zeroInnerCovariance[1],LibUtilities::ReduceMin);
        m_comm->AllReduce(zeroInnerCovariance[2],LibUtilities::ReduceMax);
        if(!zeroInnerCovariance[0])zeroInnerCovariance[1]=zeroInnerCovariance[2]=0;
        result["inner_zero_covariance"]=zeroInnerCovariance;
    }
    if(!admissible)return result;
    Vector continuity(5);for(const auto &q:Quadrature(0))
    {
        double div=0;for(int a=0;a<m_dimension;++a)div+=q.Gradient(s,offset[a])[PhysicalAxis(a)];
        const int region=q.coordinate[2]<m_parameters.geometry.height?0:1;
        continuity[region]+=q.weight*div*div;continuity[2+region]+=q.weight*div;
    }
    continuity[0]=std::sqrt(continuity[0]);continuity[1]=std::sqrt(continuity[1]);
    for(int axis:{0,1,2})if(m_dimension==3 || axis!=1)for(int side:{-1,1})
        for(const auto &f:Surface(0,axis,side<0?lower[axis]:upper[axis],side,8))
            for(int a=0;a<m_dimension;++a)continuity[4]+=f.weight*f.test.Value(s,offset[a])*f.area[PhysicalAxis(a)];
    result["strong_continuity"] = continuity;
    Vector face(8);
    for(const auto &f:Surface(0,2,m_parameters.geometry.height,1,8))
    {
        auto p=f.test.coordinate;const double measure=std::hypot(f.area[0],f.area[1],f.area[2]);
        RANSVector n{};for(int a=0;a<3;++a)n[a]=f.area[a]/measure;
        RANSTensor gi{},go{};for(int a=0;a<m_dimension;++a){gi[PhysicalAxis(a)]=Point(0,p,-1).Gradient(s,offset[a]);go[PhysicalAxis(a)]=Point(0,p,1).Gradient(s,offset[a]);}
        const double k=Point(K(),p).Value(s,offset[K()]),eps=Point(Epsilon(),p,1).Value(s,offset[Epsilon()]);
        const double pi=Point(InnerPressure(),p).Value(s,offset[InnerPressure()]),po=Point(OuterPressure(),p).Value(s,offset[OuterPressure()]);
        const auto profile=StateProfile(s,p,time,data);const double nt=profile->Original().EddyViscosity(Distance(p));
        const auto ti=RANSWallInnerMeanStress(InnerStress(k,nt,gi),m_parameters.geometry.molecularNu,pi);
        const auto c=OuterClosure(k,eps,go,p);
        RANSVector jump{};for(int a=0;a<3;++a)for(int b=0;b<3;++b)
            jump[a]+=(ti[a][b]+c.reynoldsStress[a][b]-(a==b?2*k/3-po:0))*n[b];
        double normal=0,norm=0;for(int a=0;a<3;++a){normal+=jump[a]*n[a];norm+=jump[a]*jump[a];}
        const double w=f.weight*measure;face[0]+=w;face[1]+=w*(pi-po);face[2]+=w*normal*normal;
        face[3]+=w*norm;face[4]=std::max(face[4],std::abs(normal));
        const auto gki=Point(K(),p,-1).Gradient(s,offset[K()]),gko=Point(K(),p,1).Gradient(s,offset[K()]);
        double flux=0,un=0;for(int a=0;a<3;++a){flux+=((m_parameters.geometry.molecularNu+nt/m_parameters.model.Parameters().sigmaK)*gki[a]-c.kDiffusion*gko[a])*n[a];}
        for(int a=0;a<m_dimension;++a)un+=Point(0,p).Value(s,offset[a])*n[PhysicalAxis(a)];
        face[5]+=w*flux*flux;face[6]+=w*un;face[7]+=w*normal;
    }
    face[2]=std::sqrt(face[2]/face[0]);face[3]=std::sqrt(face[3]/face[0]);face[5]=std::sqrt(face[5]/face[0]);result["interface"]=face;
    const auto raw=Residual(s,history,alpha,time,data,false,false);
    const auto constraints=Constraints(s,time,data);
    result["coupled_norms"]=Norms(Residual(s,history,alpha,time,data),constraints,ustar);
    if(NewModels())
    {
        // Partition the actual common-CG equations by basis support. A shared
        // interface row is counted once; it is not a pointwise flux jump or
        // an extra two-sided equation. No reference residual is subtracted.
        const auto physical=Residual(s,history,alpha,time,data,false);
        std::set<int> innerIds,outerIds;
        for(const auto &q:Quadrature(K()))
        {
            auto &ids=q.coordinate[2]<m_parameters.geometry.height?innerIds:outerIds;
            ids.insert(q.ids.begin(),q.ids.end());
        }
        Vector dual(3),largest(3),constantTest(3),count(3);
        for(int group=0;group<3;++group)
        {
            std::vector<int> ids;Vector load;
            for(int id=0;id<size[K()];++id)
            {
                const int kind=innerIds.count(id)&&outerIds.count(id)?2:innerIds.count(id)?0:1;
                if(kind!=group || constraints.count(offset[K()]+id))continue;
                ids.push_back(id);load.push_back(physical[offset[K()]+id]);
                largest[group]=std::max(largest[group],std::abs(load.back()));
                constantTest[group]+=m_constants[K()][id]*load.back();
            }
            count[group]=ids.size();
            if(ids.empty())continue;
            Vector inverse;
            if(SerialSparse())inverse=RANSWallAlgebra::SparseFactor(m_sparseMass[K()].Restrict(ids),true).Apply(load);
            else
            {
                WallMatrix mass(ids.size(),Vector(ids.size()));
                for(int i=0;i<ids.size();++i)for(int j=0;j<ids.size();++j)mass[i][j]=m_mass[K()][ids[i]][ids[j]];
                inverse=Factor(mass).Apply(load);
            }
            long double square=0;for(int i=0;i<load.size();++i)square+=(long double)load[i]*inverse[i];
            dual[group]=std::sqrt(std::max(0.L,square))/(K0*ustar/m_parameters.geometry.height*std::sqrt(m_volumes[K()]));
        }
        result["native_k_inner_outer_shared_interface_dual"]=dual;
        result["native_k_inner_outer_shared_interface_max_row"]=largest;
        result["native_k_inner_outer_shared_interface_constant_test"]=constantTest;
        result["native_k_inner_outer_shared_interface_row_counts"]=count;
        double velocityJump=0,kJump=0;
        for(const auto &f:Surface(0,2,m_parameters.geometry.height,1,8))
        {
            const auto &p=f.test.coordinate;
            for(int a=0;a<m_dimension;++a)velocityJump=std::max(velocityJump,
                std::abs(Point(a,p,-1).Value(s,offset[a])-Point(a,p,1).Value(s,offset[a])));
            kJump=std::max(kJump,std::abs(Point(K(),p,-1).Value(s,offset[K()])-Point(K(),p,1).Value(s,offset[K()])));
        }
        double dependent=0;for(int id:m_traceIds)dependent=std::max(dependent,std::abs(physical[offset[Epsilon()]+id]));
        result["common_velocity_k_dependent_epsilon_trace_defects"]={velocityJump,kJump,dependent};
    }
    Vector balance(2*(m_dimension+1));
    for(int b=0;b<=K();++b)
    {
        for(int i=0;i<size[b];++i){const int row=offset[b]+i;balance[b]+=m_constants[b][i]*raw[row];
            if(constraints.count(row))balance[m_dimension+1+b]+=m_constants[b][i]*raw[row];}
    }
    result["momentum_k_load_and_boundary_reaction"]=balance;
    Vector wall(5),repair(8);repair[6]=1;
    const bool saturated=m_innerStress==RANSWallInnerStressModel::SaturatedR2;
    // R2-only layouts: [active physical volume, positive-k volume, min a/c,
    // max a/c, zero-k physical volume, minimum turbulent radial tangent].
    Vector innerLaw{0,0,sentinel,-sentinel,0,sentinel},innerChange(3),ratioPoint(10);
    double localMaximumRatio=-sentinel;
    for(const auto &ray:m_rays)for(const auto &q:ray.points)
    {
        RANSTensor G{};for(int a=0;a<m_dimension;++a)G[PhysicalAxis(a)]=q.v.Gradient(s,offset[a]);
        const double k=q.v.Value(s,offset[K()]),nt=ray.profile->Original().EddyViscosity(Distance(q.v.coordinate));
        const auto selected=InnerStress(k,nt,G);const double production=selected.production;
        double strain=0;for(const auto &row:selected.strain)for(double value:row)strain+=value*value;
        if(saturated)
        {
            const auto linear=EvaluateRANSWallInnerStress(k,nt,G);
            for(int a=0;a<3;++a)for(int b=0;b<3;++b)
                innerChange[0]-=q.weight*(selected.reynoldsStress[a][b]-linear.reynoldsStress[a][b])*G[a][b];
            innerChange[1]+=q.weight*(selected.production-linear.production);
            if(selected.active)innerLaw[0]+=q.weight;
            if(k>0)
            {
                const double ratio=selected.a/selected.c;innerLaw[1]+=q.weight;
                innerLaw[2]=std::min(innerLaw[2],ratio);innerLaw[3]=std::max(innerLaw[3],ratio);
                if(ratio>localMaximumRatio)
                {localMaximumRatio=ratio;ratioPoint={q.v.coordinate[0],q.v.coordinate[1],q.v.coordinate[2],
                    k,ray.profile->K(),nt,selected.strainNorm,ratio,selected.t,selected.radialTangent};}
            }
            else innerLaw[4]+=q.weight;
            innerLaw[5]=std::min(innerLaw[5],selected.radialTangent);
        }
        const auto pre=ray.profile->AtProduction(Distance(q.v.coordinate),production);
        const auto values=WallSource(s,q.v.coordinate,production,*ray.profile,m_cappedDestruction?&ray.destruction:nullptr);
        wall[0]+=q.weight*ray.profile->Original().OriginalDissipation(Distance(q.v.coordinate));wall[1]+=q.weight*values.dissipation;
        wall[2]+=q.weight*values.redistribution;wall[3]+=q.weight*production;
        wall[4]+=q.weight*m_parameters.geometry.molecularNu*strain/2;
        repair[0]+=q.weight*pre.dissipation;repair[1]+=q.weight*values.dissipation;
        repair[2]+=q.weight*(values.dissipation-values.redistribution);repair[3]+=q.weight*values.redistribution;
        if(m_cappedDestruction)
        {
            const auto point=EvaluateRANSCappedDestructionPoint(k,ray.profile->K(),pre.dissipation,
                pre.dissipation-pre.redistribution,ray.destruction.b,ray.destruction.m,m_destruction);
            if(point.limiter<1)repair[4]+=q.weight;
            if(ray.destruction.capped)repair[5]+=q.weight;
            repair[6]=std::min(repair[6],ray.destruction.lambda);
            repair[7]=std::max(repair[7],point.rate);
        }
    }
    m_comm->AllReduce(wall,LibUtilities::ReduceSum);
    if(AcceptedAWF())result["accepted_awf_sources"]={wall[3],wall[1],wall[2],wall[4],wall[3]-wall[1]};
    else result["wall_sources"]=wall;
    for(int i=0;i<6;++i)m_comm->AllReduce(repair[i],LibUtilities::ReduceSum);
    m_comm->AllReduce(repair[6],LibUtilities::ReduceMin);m_comm->AllReduce(repair[7],LibUtilities::ReduceMax);
    if(!AcceptedAWF())result["robust_source"]=repair;
    if(saturated)
    {
        for(int i:{0,1,4})m_comm->AllReduce(innerLaw[i],LibUtilities::ReduceSum);
        for(int i:{2,5})m_comm->AllReduce(innerLaw[i],LibUtilities::ReduceMin);
        m_comm->AllReduce(innerLaw[3],LibUtilities::ReduceMax);
        int rank=innerLaw[1]>0 && localMaximumRatio==innerLaw[3]?m_comm->GetRank():m_comm->GetSize();
        m_comm->AllReduce(rank,LibUtilities::ReduceMin);
        if(rank<m_comm->GetSize())m_comm->Bcast(ratioPoint,rank);
        for(double &value:innerLaw)if(value==sentinel || value==-sentinel)value=0;
        m_comm->AllReduce(innerChange,LibUtilities::ReduceSum);innerChange[2]=innerChange[0]-innerChange[1];
        result["inner_stress"]=innerLaw;result["inner_stress_max_ratio_point"]=ratioPoint;
        result["inner_constitutive_modification"]=innerChange;
    }
    Vector modification(3);
    if(m_variantName!="Original")for(const auto &q:m_outer)
    {
        RANSTensor G{};for(int a=0;a<m_dimension;++a)G[PhysicalAxis(a)]=q.v.Gradient(s,offset[a]);
        const double k=q.v.Value(s,offset[K()]),eps=q.e.Value(s,offset[Epsilon()]);
        const auto original=EvaluateRANSClosure(m_parameters.model.Model(),k,eps,G,Distance(q.v.coordinate),m_parameters.model.Parameters(),0);
        const auto changed=OuterClosure(k,eps,G,q.v.coordinate);
        for(int a=0;a<3;++a)for(int b=0;b<3;++b)
            modification[0]-=q.weight*(changed.reynoldsStress[a][b]-original.reynoldsStress[a][b])*G[a][b];
        modification[1]+=q.weight*(changed.production-original.production);
    }
    m_comm->AllReduce(modification,LibUtilities::ReduceSum);
    // Mean weak stress work equals production, and the k residual carries
    // minus production. Their sum in total stored energy is zero.
    modification[2]=modification[0]-modification[1];result["constitutive_modification"]=modification;
    auto energy=EnergyBudget(s,history,alpha,time,data);double reaction=0;
    for(int b=0;b<=K();++b)for(int i=0;i<size[b];++i)if(constraints.count(offset[b]+i))
        reaction+=(b==K()?m_constants[b][i]:s[offset[b]+i])*raw[offset[b]+i];
    energy.push_back(reaction);energy.push_back(AnchorEnergyWork(s));energy.push_back(energy[0]-reaction-energy[10]);
    result["energy"]=energy;result["roughness"]=RoughnessBudget(s,time,data);
    return result;
}
std::string RANSWallEvolutionPatch::InitializationReceipt(const Vector &state,const Vector &reference,double time,
    const RANSWallEvolutionData &data,const RANSWallEvolutionData &implicit,const std::map<int,double> &generated) const
{
    namespace H=RANSWallR3;namespace J=RANSWallSolveJson;
    const auto live=ExternalConstraints(time,data),other=ExternalConstraints(time,implicit);
    auto hash=[](const std::map<int,double> &map){Vector bytes;for(const auto &[row,value]:map){bytes.push_back(row);bytes.push_back(value);}return H::Hash(bytes);};
    const auto effective=ExternalLift(state,time,data);double difference=0,boundaryDifference=0;
    for(int i=0;i<state.size();++i)difference=std::max(difference,std::abs(state[i]-effective[i]));
    for(const auto &[row,value]:live)if(generated.count(row))boundaryDifference=std::max(boundaryDifference,std::abs(value-generated.at(row)));
    std::vector<std::pair<std::string,std::string>> entries{{"raw_initial_sha256",J::Quote(H::Hash(state))},{"effective_initial_sha256",J::Quote(H::Hash(effective))},
        {"fixed_reference_sha256",J::Quote(H::Hash(reference))},{"explicit_boundary_sha256",J::Quote(hash(live))},
        {"implicit_boundary_sha256",J::Quote(hash(other))},{"boundary_hash_format",J::Quote("sorted row-as-f64,value pairs, native little endian")},
        {"raw_effective_max_difference",J::Number(difference)},{"source_target_boundary_max_difference",J::Number(boundaryDifference)},
        {"live_fixed_boundary",J::Constraints(live)},{"effective_initial",J::Numbers(effective)}};
    if(OpenBoundary())
    {
        const auto explicitLoads=OpenBoundaryContributions(state,time,data),implicitLoads=OpenBoundaryContributions(state,time,implicit);
        Require(data.openBoundaryIdentity==implicit.openBoundaryIdentity && explicitLoads.traction==implicitLoads.traction &&
            explicitLoads.diffusion==implicitLoads.diffusion,"explicit/implicit open boundary data differ");
        entries.push_back({"open_boundary_definition",J::Quote(data.openBoundaryIdentity)});
        entries.push_back({"open_traction_load_sha256",J::Quote(H::Hash(explicitLoads.traction))});
        entries.push_back({"open_diffusive_load_sha256",J::Quote(H::Hash(explicitLoads.diffusion))});
        entries.push_back({"open_advective_load_sha256",J::Quote(H::Hash(explicitLoads.advection))});
        entries.push_back({"open_explicit_implicit_agreement","true"});
    }
    return J::Object(entries);
}
std::string RANSWallEvolutionPatch::PreservationIdentity(bool legacyOperator) const
{
    std::ostringstream out;out<<std::setprecision(17)<<"NeutralVerification-v1:"<<m_parameters.model.Signature()
        <<":"<<m_parameters.geometry.height<<":"<<m_parameters.geometry.roughness<<":"<<m_parameters.geometry.sampleHeight
        <<":"<<m_parameters.geometry.tailStart<<":"<<m_parameters.chi<<":"<<WallPolicyName()<<":"<<m_robustUstar;
    for(int b=0;b<Gauge();++b)
    {
        out<<"|"<<offset[b]<<","<<size[b];
        for(const auto &e:fields[b]->Elements())out<<":"<<e.geometry<<","<<e.degree<<","<<e.lower[0]<<","<<e.upper[0]<<","<<e.lower[2]<<","<<e.upper[2];
        for(const auto &q:fields[b]->Quadrature(2)){for(int id:q.ids)out<<","<<id;}
        if(fields[b]->Elements().front().Directional(m_dimension))
        {
            Require(!legacyOperator,"directional layout has no legacy ownership identity");
            out<<":UniformDirectionalModifiedA2DV1";
            // Deliberately exclude points/quadrature: frozen coefficients may
            // be evaluated with independently increased integration orders.
            for(int e=0;e<fields[b]->Field()->GetExpSize();++e)
                for(int d=0;d<m_dimension;++d)
                    out<<':'<<fields[b]->Field()->GetExp(e)->GetBasis(d)->GetBasisType()
                       <<':'<<fields[b]->Field()->GetExp(e)->GetBasisNumModes(d);
        }
    }
    if(!legacyOperator)out<<":numerical_operator="<<RANSWallNumericalOperatorIdentity();
    if(SerialSparse())out<<":linear_backend=AcceptedSerialSparseV1:CSR1:EigenLU-LDLT:FGMRES80x1000:relative2e-12:absolute1e-25:consistent-free-M:lumped-pressure-preconditioner-only";
    if(CurvedVerification())out<<":curved="<<CurvedCapabilityName()<<":"<<m_mapping.Signature()<<":CT6-conservative-v1";
    if(OpenBoundary())
    {
        Require(!legacyOperator,"open boundary has no legacy replay identity");
        out<<":boundary="<<m_openBoundaryIdentity;
    }
    return out.str();
}
std::string RANSWallEvolutionPatch::PreservationPacket(const Vector &state,double time,
    const RANSWallEvolutionData &data,bool basis,const std::vector<RANSWallCoordinate> &grid,const std::vector<int> &sides,bool ct6Primitives) const
{
    namespace J=RANSWallSolveJson;
    Require(NewModels() && Dimension()==2 && (!Mapped() || CurvedVerification()) && m_comm->GetSize()==1,
        "verification is model-aware serial flat2D only");
    Require(m_variantName=="Original" && !m_cappedDestruction && m_innerStress==RANSWallInnerStressModel::Linear,
        "verification requires original outer, linear inner and prescribed uncapped sources");
    EvaluationScope evaluation(*this,data,time);Refresh(state,time,&data);
    auto row=[&](const RANSWallGlobalPoint &p,int block=-1)
    {std::vector<std::pair<std::string,std::string>> result{{"ids",J::Numbers(p.ids)},{"phi",J::Numbers(p.phi)},
        {"constant",J::Numbers(p.constant)},
        {"derivative",J::Rows(p.derivative)},{"geometry",J::Number(p.geometry)},
        {"element",J::Number(p.element)},
        {"endpoint_roundoff_axes",J::Number(p.endpointRoundoffAxes)},
        {"raw_local_coordinate",J::Numbers(p.rawLocalCoordinate)}};
     if(CurvedVerification() && ct6Primitives && block>=0)
     {
         const auto chart=fields[block]->ElementPoint(p.element,p.coordinate);
         Require(chart.ids==p.ids,"CT6 exported chart/physical row ownership differs");
         result.push_back({"native_chart_derivative",J::Rows(chart.derivative)});
         result.push_back({"native_chart_phi",J::Numbers(chart.phi)});
     }
     return J::Object(result);};
    std::unique_ptr<RANSWallManufactured> manufactured;
    if(CurvedManufactured() && ct6Primitives)manufactured=std::make_unique<RANSWallManufactured>(*this);
    auto point=[&](const VolumePoint &q,int region)
    {
        std::vector<std::pair<std::string,std::string>> fields{{"coordinate",J::Numbers(q.v.coordinate)},{"weight",J::Number(q.weight)},
            {"region",J::Number(region)},{"velocity_k",row(q.v,0)},{"pressure",row(q.p,InnerPressure()+region)},
            {"physical_coordinate",J::Numbers(Physical(q.v.coordinate))},
            {"volume_metric",J::Number(m_mapping.At(q.v.coordinate).determinant)},
            {"normal_J",J::Number(m_mapping.At(q.v.coordinate).J)},
            {"generating_geometry",J::Number(q.generatingGeometry)},
            {"generating_radial_panel",J::Numbers(Vector{q.panelLower,q.panelUpper})},
            {"epsilon",region==1?row(q.e,Epsilon()):"null"}};
        if(manufactured)
        {
            const auto &chart=q.v.coordinate;auto H=chart;H[2]=m_parameters.geometry.height;
            const auto metric=m_mapping.At(H);
            const auto exact=manufactured->At(chart,time,region==1);
            Vector body,outerBody;
            for(int a=0;a<m_dimension;++a)
            {
                body.push_back((region?data.outerBody:data.innerBody)[a](chart,time));
                outerBody.push_back(data.outerBody[a](H,time));
            }
            fields.push_back({"CT6_prescribed_data",J::Object({
                {"body",J::Numbers(body)},{"outer_body_at_H",J::Numbers(outerBody)},
                {"k_source",J::Number((region?data.outerKSource:data.innerKSource)(chart,time))},
                {"epsilon_source",region?J::Number(data.outerEpsilonSource(chart,time)):"null"},
                {"H_velocity_k_row",row(Point(0,H,1),0)},
                {"H_pressure_row",row(Point(OuterPressure(),H,1),OuterPressure())},
                {"H_epsilon_row",row(Point(Epsilon(),H,1),Epsilon())},
                {"H_normal",J::Numbers(metric.normal)},{"H_Br",J::Rows(metric.Br)},
                {"continuous_values",J::Numbers(exact.value)},
                {"continuous_gradients",J::Rows(exact.gradient)},
                {"continuous_R",J::Rows(exact.reynolds)},{"continuous_Q",J::Rows(exact.quadratic)},
                {"continuous_T",J::Rows(exact.stress)},
                {"continuous_dtU",J::Numbers(exact.velocityTime)},
                {"continuous_dtk",J::Number(exact.kTime)},
                {"continuous_dtepsilon",J::Number(exact.epsilonTime)},
                {"continuous_required_body",J::Numbers(exact.body)},
                {"continuous_Dk",J::Number(exact.kDiffusion)},
                {"continuous_Depsilon",J::Number(exact.epsilonDiffusion)},
                {"continuous_production",J::Number(exact.production)}})});
        }
        return J::Object(fields);
    };
    std::vector<std::string> volume,bodyVolume,faces,cellFaces,mass,profiles,values,gradients,regions,observedSides,provenance,tractions,dependentEpsilon;
    std::size_t ni=0,no=m_outer.size();
    for(const auto &ray:m_rays)
    {
        ni+=ray.points.size();const auto &ref=ray.profile->Reference();
        profiles.push_back(J::Object({{"coordinate",J::Numbers(ray.coordinate)},{"K",J::Number(ray.profile->K())},
            {"chi",J::Number(ref.weight)},{"cutoff",J::Number(ref.cutoff)},{"budget_ratio",J::Number(ref.productionBudget/ref.minimumBudget)},
            {"viscosity_intercept",J::Number(ray.profile->Original().ViscosityInterceptHeight())},
            {"source_panels",J::Numbers(m_options.physicalBreakpoints)}}));
        if(basis)for(const auto &q:ray.points)volume.push_back(point(q,0));
    }
    if(basis)for(const auto &q:m_outer)volume.push_back(point(q,1));
    if(basis && manufactured)
    {
        // Prescribed IMEX body loads retain their own fixed generating panels;
        // do not substitute the current-K wall-source quadrature for them.
        for(const auto &q:m_bodyInner)bodyVolume.push_back(point(q,0));
        for(const auto &q:m_outer)bodyVolume.push_back(point(q,1));
    }
    const double H=m_parameters.geometry.height;
    if(basis)for(int region=0;region<2;++region)
    {
        for(int axis:{0,2})for(int side:{-1,1})
        {
            const double location=axis==0?(side<0?lower[0]:upper[0]):
                (side<0?(region==0?0.:H):(region==0?H:upper[2]));
            for(const auto &f:Surface(0,axis,location,side,m_interfaceOrder?m_interfaceOrder:8))
            {
                if((f.test.coordinate[2]<H?0:1)!=region && axis==0)continue;
                if(axis==2 && location==H)
                {
                    const auto &e=fields[0]->Elements().at(f.test.element);
                    if((e.lower[2]>=H?1:0)!=region)continue;
                }
                const auto pressure=GeometryPoint(InnerPressure()+region,f.test.geometry,f.test.coordinate);
                std::vector<std::pair<std::string,std::string>> record{{"coordinate",J::Numbers(f.test.coordinate)},{"weight",J::Number(f.weight)},
                    {"area",J::Numbers(f.area)},{"region",J::Number(region)},
                    {"interface",J::Boolean(axis==2 && location==H)},
                    {"velocity_k",row(f.test,0)},{"pressure",row(pressure,InnerPressure()+region)},
                    {"epsilon",region==1?row(GeometryPoint(Epsilon(),f.test.geometry,f.test.coordinate),Epsilon()):"null"}};
                if(manufactured)
                {
                    const auto exact=manufactured->At(f.test.coordinate,time,region==1);
                    record.push_back({"CT6_continuous",J::Object({{"values",J::Numbers(exact.value)},
                        {"gradients",J::Rows(exact.gradient)},{"T",J::Rows(exact.stress)},
                        {"Dk",J::Number(exact.kDiffusion)},{"Depsilon",J::Number(exact.epsilonDiffusion)}})});
                }
                faces.push_back(J::Object(record));
            }
        }
        for(const auto &q:Quadrature(InnerPressure()+region))mass.push_back(J::Object({
            {"coordinate",J::Numbers(q.coordinate)},{"weight",J::Number(q.weight)},{"region",J::Number(region)},{"pressure",row(q,InnerPressure()+region)}}));
    }
    if(basis && CurvedVerification() && ct6Primitives)for(int block:{0,Epsilon(),InnerPressure(),OuterPressure()})
    {
        std::array<std::set<double>,3> locations;
        for(const auto &e:fields[block]->Elements())for(int axis:{0,2})
        {locations[axis].insert(e.lower[axis]);locations[axis].insert(e.upper[axis]);}
        for(int axis:{0,2})for(double location:locations[axis])for(int side:{-1,1})
            for(const auto &f:Surface(block,axis,location,side,m_interfaceOrder?m_interfaceOrder:12))
                cellFaces.push_back(J::Object({{"block",J::Number(block)},{"axis",J::Number(axis)},
                    {"side",J::Number(side)},{"coordinate",J::Numbers(f.test.coordinate)},
                    {"weight",J::Number(f.weight)},{"area",J::Numbers(f.area)},
                    {"row",row(f.test,block)}}));
    }
    Require(sides.empty() || sides.size()==grid.size(),"sample side count differs");
    for(int index=0;index<grid.size();++index)
    {
        const auto &p=grid[index];const bool interface=p[2]==H;
        const int side=sides.empty()?0:sides[index];
        Require(interface?(side==-1 || side==1):side==0,"ambiguous interface sample: declare inner=-1 or outer=+1 at H; interior side must be zero");
        const int region=interface?(side<0?0:1):(p[2]<H?0:1);Vector v,geometry;std::vector<Vector> g;
        for(int b:{0,1,K(),Epsilon(),InnerPressure()+region})
        {
            if(b==Epsilon() && region==0){v.push_back(0.);g.push_back({0,0,0});geometry.push_back(-1);continue;}
            const auto q=Point(b,p,side);v.push_back(q.Value(state,offset[b]));geometry.push_back(q.geometry);
            const auto derivative=q.Gradient(state,offset[b]);g.emplace_back(derivative.begin(),derivative.end());
        }
        RANSTensor G{},stress{};for(int a=0;a<2;++a)for(int j=0;j<3;++j)G[PhysicalAxis(a)][j]=g[a][j];
        if(region==0){const auto profile=StateProfile(state,p,time,data);stress=RANSWallInnerMeanStress(InnerStress(v[2],profile->Original().EddyViscosity(Distance(p)),G),m_parameters.geometry.molecularNu,v[4]);}
        else{const auto law=OuterClosure(v[2],v[3],G,p);for(int a=0;a<3;++a)for(int j=0;j<3;++j)stress[a][j]=-law.reynoldsStress[a][j]+(a==j?2*v[2]/3-v[4]:0.);}
        Vector traction;auto normal=m_mapping.At(p).normal;
        if(CurvedVerification() && interface){normal=m_mapping.MatchingArea(p);const double norm=std::hypot(normal[0],normal[1],normal[2]);for(double &v:normal)v/=norm;}
        for(int a=0;a<3;++a){double v=stress[a][2];if(CurvedVerification()){v=0;for(int b=0;b<3;++b)v+=stress[a][b]*normal[b];}traction.push_back(v);}tractions.push_back(J::Numbers(traction));
        dependentEpsilon.push_back(region==0?J::Number(DependentEpsilon(state,p,time,data)):"null");
        values.push_back(J::Numbers(v));gradients.push_back(J::Rows(g));regions.push_back(J::Number(region));
        observedSides.push_back(J::Number(side));provenance.push_back(J::Numbers(geometry));
    }
    const auto audit=CompleteRestrictedAudit(state,state,0,time,data);
    std::vector<std::pair<std::string,std::string>> diagnostic;
    for(const auto &[key,v]:audit)diagnostic.push_back({key,J::Numbers(v)});
    std::vector<std::string> elements;
    for(const auto &e:fields[0]->Elements())elements.push_back(J::Object({{"geometry",J::Number(e.geometry)},
        {"degree",J::Number(e.degree)},{"degrees",J::Numbers(e.degrees)},
        {"lower",J::Numbers(e.lower)},{"upper",J::Numbers(e.upper)}}));
    std::vector<std::string> fieldBasis;
    for(int b=0;b<Gauge();++b)
    {
        std::vector<std::string> nativeElements;
        for(int e=0;e<fields[b]->Field()->GetExpSize();++e)
        {
            const auto exp=fields[b]->Field()->GetExp(e);std::vector<std::string> axes;
            for(int d=0;d<m_dimension;++d)
                axes.push_back(J::Numbers(Vector{double(exp->GetBasis(d)->GetBasisType()),
                    double(exp->GetBasisNumModes(d)),double(exp->GetBasis(d)->GetPointsType()),
                    double(exp->GetBasis(d)->GetNumPoints())}));
            nativeElements.push_back(J::Object({{"geometry",J::Number(exp->GetGeom()->GetGlobalID())},
                {"axes_type_modes_pointsType_points",J::List(axes)}}));
        }
        fieldBasis.push_back(J::List(nativeElements));
    }
    return J::Object({{"schema","2"},{"read_only","true"},{"time",J::Number(time)},{"identity",J::Quote(PreservationIdentity())},
        {"numerical_operator",J::Quote(RANSWallNumericalOperatorIdentity())},
        {"resolved_descriptor",J::Quote(m_parameters.model.Signature())},{"state",J::Numbers(state)},
        {"offsets",J::Numbers(offset)},{"sizes",J::Numbers(size)},{"elements",J::List(elements)},
        {"field_basis_configuration",J::List(fieldBasis)},
        {"actual_quadrature",J::Object({{"wall_normal",J::Number(m_parameters.normalOrder)},
            {"wall_tangential",J::Number(m_parameters.tangentialOrder)},
            {"outer_volume",J::Number(m_volumeOrder?m_volumeOrder:m_parameters.tangentialOrder)},
            {"mass_override",J::Number(m_massOrder)},{"mass_default",J::Quote("max(8,p+3) per native field")},
            {"interface",J::Number(m_interfaceOrder?m_interfaceOrder:8)},
            {"reference_projection",J::Quote(OpenBoundary()?
                "Functions: native per-field GLL FwdTrans seed, then three consistent-mass residual corrections with default Gauss max(8,p+3) per axis for nonconstant data":
                "Gauss32 on native radial elements, supplied breakpoints and geometric 1e-4*1.5^j subpanels")},
            {"inner_volume_points",J::Number(ni)},{"outer_volume_points",J::Number(no)}})},
        {"profiles",J::List(profiles)},{"volume",J::List(volume)},{"CT6_body_volume",J::List(bodyVolume)},
        {"CT6_cell_faces",J::List(cellFaces)},
        {"faces",J::List(faces)},{"pressure_mass",J::List(mass)},
        {"audit_alpha_zero",J::Object(diagnostic)},
        {"grid",J::Rows(grid)},{"grid_regions",J::List(regions)},
        {"grid_inner_prescribed_epsilon",J::List(dependentEpsilon)},{"grid_sides",J::List(observedSides)},{"grid_geometry_by_field",J::List(provenance)},{"grid_traction",J::List(tractions)},{"grid_values",J::List(values)},{"grid_gradients",J::List(gradients)}});
}

namespace
{
namespace CJ=RANSWallSolveJson;
using CV=RANSWallAlgebra::Vector;
struct Spectrum { CV eigenvalues; RANSTensor vectors{}; };
Spectrum CovarianceSpectrum(const RANSTensor &a)
{
    const auto eig=CJ::Eigensystem(a);
    return {CV(eig.eigenvalues.begin(),eig.eigenvalues.end()),eig.vectors};
}
std::string TensorJSON(const RANSTensor &a)
{std::vector<std::string> rows;for(const auto &r:a)rows.push_back(CJ::Numbers(r));return CJ::List(rows);}
std::string EvaluationRowJSON(const RANSWallGlobalPoint &p)
{return CJ::Object({{"geometry",CJ::Number(p.geometry)},{"element",CJ::Number(p.element)},
    {"ids",CJ::Numbers(p.ids)},{"phi",CJ::Numbers(p.phi)},{"constant",CJ::Numbers(p.constant)},
    {"derivative",CJ::List({CJ::Numbers(p.derivative[0]),CJ::Numbers(p.derivative[1]),CJ::Numbers(p.derivative[2])})}});}
std::string LawJSON(const RANSClosureValues &c,double k,const RANSTensor &g,const std::string &remainder="regularization")
{
    const bool inner=remainder=="inner_stress_remainder";
    const auto eig=CovarianceSpectrum(c.reynoldsStress);RANSVector direction{};
    for(int a=0;a<3;++a)direction[a]=eig.vectors[a][0];
    RANSTensor iso{},linear{},reg{};
    for(int a=0;a<3;++a)for(int b=0;b<3;++b)
    {iso[a][b]=a==b?2*k/3:0;linear[a][b]=-c.eddyViscosity*(g[a][b]+g[b][a]);
     reg[a][b]=c.reynoldsStress[a][b]-iso[a][b]-linear[a][b]-c.quadraticStress[a][b]-c.cubicStress[a][b];}
    std::vector<std::pair<std::string,std::string>> terms,projections;
    for(const auto &[name,tensor]:std::map<std::string,RANSTensor>{{"isotropic",iso},{"linear",linear},{"quadratic",c.quadraticStress},{"cubic",c.cubicStress},{remainder,reg}})
    {double value=0;for(int a=0;a<3;++a)for(int b=0;b<3;++b)value+=direction[a]*tensor[a][b]*direction[b];terms.push_back({name,TensorJSON(tensor)});projections.push_back({name,CJ::Number(value)});}
    return CJ::Object({{"tensor",TensorJSON(c.reynoldsStress)},{"eigenvalues",CJ::Numbers(eig.eigenvalues)},
        {"critical_eigenvector",CJ::Numbers(direction)},{"terms",CJ::Object(terms)},
        {"critical_term_contractions",CJ::Object(projections)},{"c_mu",CJ::Number(c.cMu)},
        {"eddy_viscosity",CJ::Number(c.eddyViscosity)},{"strain_invariant",CJ::Number(c.strainInvariant)},
        {"rotation_invariant",CJ::Number(c.rotationInvariant)},{"invariant_branch",CJ::Quote(inner?"not_applicable_to_inner_law":(c.strainInvariant>=c.rotationInvariant?"strain":"rotation"))},
        {"production",CJ::Number(c.production)},{"k_source",CJ::Number(c.kSource)},
        {"epsilon_production",inner?"null":CJ::Number(c.epsilonProduction)},{"epsilon_source",inner?"null":CJ::Number(c.epsilonSource)},
        {"yap",inner?"null":CJ::Number(c.yap)},{"k_diffusion",CJ::Number(c.kDiffusion)},{"epsilon_diffusion",inner?"null":CJ::Number(c.epsilonDiffusion)}});
}
}
std::string RANSWallEvolutionPatch::ConstitutiveStateDiagnostic(const Vector &s,double time,
    const RANSWallEvolutionData &data,const std::vector<RANSWallCoordinate> &witnesses) const
{
    Require(!NewModels(),"CLS-specific diagnostic/action unsupported for A93/A95/E94E95; use native model-aware audits");
    namespace J=RANSWallSolveJson;Require(AcceptedAWF() && !Mapped() && m_dimension==2,"constitutive forensic scope is accepted AWF flat2D");
    const auto before=RANSWallR3::Hash(s);const auto domain=ScalarDomain(s,time,data,true);
    if(!domain.admissible)return J::Object({{"scalar_budget",ScalarDomainJSON(domain)}, {"covariance_evaluated","false"},{"state_sha256",J::Quote(before)}});
    Refresh(s,time,&data);double minimum[4]={1e300,1e300,1e300,1e300};std::string minima[4]={"null","null","null","null"};std::map<std::string,int> counts;
    auto evaluate=[&](const RANSWallGlobalPoint &v,bool inner,const std::string &family,bool retain,double epsilonOverride=std::numeric_limits<double>::quiet_NaN())
    {
        RANSTensor g{};for(int a=0;a<m_dimension;++a)g[PhysicalAxis(a)]=v.Gradient(s,offset[a]);
        const double k=GeometryPoint(K(),v.geometry,v.coordinate).Value(s,offset[K()]);const auto profile=StateProfile(s,v.coordinate,time,data);
        RANSWallGlobalPoint e;
        if(!inner)
        {
            e=GeometryPoint(Epsilon(),v.geometry,v.coordinate);
        }
        const double retainedEpsilon=inner?profile->Dissipation(Distance(v.coordinate)):e.Value(s,offset[Epsilon()]);
        const double epsilon=std::isfinite(epsilonOverride)?epsilonOverride:retainedEpsilon;
        RANSClosureValues c{},original{},regularized{};
        if(inner)
        {
            const double nt=profile->Original().EddyViscosity(Distance(v.coordinate));const auto wall=InnerStress(k,nt,g);
            c.reynoldsStress=wall.reynoldsStress;c.production=wall.production;c.eddyViscosity=nt;
            c.cMu=profile->Original().Parameters().cMu;double ss=0,oo=0;
            for(int a=0;a<3;++a)for(int b=0;b<3;++b){ss+=std::pow(g[a][b]+g[b][a],2);oo+=std::pow(g[a][b]-g[b][a],2);}
            c.strainInvariant=k/epsilon*std::sqrt(ss/2);c.rotationInvariant=k/epsilon*std::sqrt(oo/2);
            for(int a=0;a<3;++a)for(int b=0;b<3;++b)c.divergence+=a==b?g[a][b]:0;
            c.kSource=c.production-epsilon;c.kDiffusion=m_parameters.geometry.molecularNu+nt/m_parameters.model.Parameters().sigmaK;
        }
        else
        {
            c=OuterClosure(k,epsilon,g,v.coordinate);original=EvaluateRANSClosure(m_parameters.model.Model(),k,epsilon,g,Distance(v.coordinate),m_parameters.model.Parameters(),0);
            RANSRobustVariant variant;variant.beta=.4;variant.gamma=m_parameters.model.Kappa()*m_robustUstar*Distance(v.coordinate);
            const double shear=m_robustUstar/(m_parameters.model.Kappa()*Distance(v.coordinate));
            for(int a=0;a<3;++a)for(int b=0;b<3;++b)variant.referenceStrain[a][b]=shear*(m_referenceFlow[a]*m_referenceNormal[b]+m_referenceNormal[a]*m_referenceFlow[b]);
            regularized=EvaluateRANSRobustClosure(m_parameters.model.Model(),k,epsilon,g,Distance(v.coordinate),m_parameters.model.Parameters(),0,variant);
        }
        const double eigen=CovarianceSpectrum(c.reynoldsStress).eigenvalues[0];const int region=inner?0:1;++counts[family];
        const double originalEigen=inner?1e300:CovarianceSpectrum(original.reynoldsStress).eigenvalues[0];
        const double regularizedEigen=inner?1e300:CovarianceSpectrum(regularized.reynoldsStress).eigenvalues[0];
        if(!retain && eigen>=minimum[region] && originalEigen>=minimum[2] && regularizedEigen>=minimum[3])return std::string("null");
        const auto kr=SampleKRow(v.coordinate);const auto epsH=Point(Epsilon(),{v.coordinate[0],0,m_parameters.geometry.height},1);
        const bool rowInner=v.element>=0 && fields[0]->Elements().at(v.element).upper[2]<=m_parameters.geometry.height+1e-12;
        const auto item=J::Object({{"coordinate",J::Numbers(v.coordinate)},{"geometry",J::Number(v.geometry)},
            {"element",J::Number(v.element)},{"inner",J::Boolean(inner)},{"inventory",J::Quote(family)},
            {"gradient_element_region",J::Quote(rowInner?"inner":"outer")},
            {"constitutive_region",J::Quote(inner?"inner":"outer")},
            {"gradient_region_matches_law",J::Boolean(rowInner==inner)},
            {"k",J::Number(k)},{"epsilon",J::Number(epsilon)},{"epsilon_retained",J::Number(retainedEpsilon)},
            {"epsilon_semantics",J::Quote(std::isfinite(epsilonOverride)?"refreshed dependent-target audit row":"actual retained field; inner epsilon dependent")},{"endpoint_K",J::Number(profile->K())},
            {"cutoff",J::Number(profile->Reference().cutoff)},{"budget_ratio",J::Number(profile->Reference().productionBudget/profile->Reference().minimumBudget)},
            {"velocity_gradient",TensorJSON(g)},{"divergence",J::Number(c.divergence)},
            {"velocity_evaluation_row",EvaluationRowJSON(v)},{"epsilon_evaluation_row",inner?"null":EvaluationRowJSON(e)},
            {"endpoint_evaluation_row",EvaluationRowJSON(kr)},{"epsilon_H_retained",J::Number(epsH.Value(s,offset[Epsilon()]))},
            {"epsilon_H_refreshed_projected",domain.traceState.empty()?"null":J::Number(epsH.Value(domain.traceState,offset[Epsilon()]))},
            {"epsilon_H_current_unprojected",J::Number(profile->Dissipation(m_parameters.geometry.height))},
            {"selected",LawJSON(c,k,g,inner?"inner_stress_remainder":"regularization")},{"same_state_Original",inner?"null":LawJSON(original,k,g)},
            {"same_state_Regularized",inner?"null":LawJSON(regularized,k,g)}});
        if(eigen<minimum[region]){minimum[region]=eigen;minima[region]=item;}
        if(originalEigen<minimum[2]){minimum[2]=originalEigen;minima[2]=item;}
        if(regularizedEigen<minimum[3]){minimum[3]=regularizedEigen;minima[3]=item;}return item;
    };
    for(const auto &point:domain.observations)
    {
        if(point.block!=K() && !point.dependentTrace)continue;
        const auto provenance=ResolveObservation(point);
        evaluate(*provenance.velocity,provenance.inner,point.inventory,false,
            point.dependentTrace?point.value:std::numeric_limits<double>::quiet_NaN());
    }
    const auto nativeInner=minima[0],nativeOuter=minima[1];
    std::vector<std::string> neighbors;
    for(const auto &point:witnesses)for(int j=0;j<fields[0]->Elements().size();++j)
    {
        const auto &element=fields[0]->Elements()[j];bool contains=true;
        for(int a:{0,2})contains=contains && point[a]>=element.lower[a]-1e-12 && point[a]<=element.upper[a]+1e-12;
        if(contains)neighbors.push_back(evaluate(fields[0]->ElementPoint(j,point),element.upper[2]<=m_parameters.geometry.height+1e-12,"explicit_one_sided_witness",true));
    }
    Require(before==RANSWallR3::Hash(s),"forensic tensor diagnostic changed selected coefficients");
    std::vector<std::pair<std::string,std::string>> inventory;for(auto [name,n]:counts)inventory.push_back({name,J::Number(n)});
    return J::Object({{"audit_revision",J::Quote(AuditRevision())},{"state_sha256",J::Quote(before)},{"scalar_budget",ScalarDomainJSON(domain)},
        {"minimum_inner",nativeInner},{"minimum_outer",nativeOuter},
        {"minimum_with_explicit_neighbors_inner",minima[0]},{"minimum_with_explicit_neighbors_outer",minima[1]},
        {"minimum_same_state_Original_outer",minima[2]},{"minimum_same_state_Regularized_outer",minima[3]},{"one_sided_witnesses",J::List(neighbors)},
        {"inventory_counts",J::Object(inventory)},{"offsets",J::Numbers(offset)},
        {"same_state_Regularized_definition",J::Quote("beta=.4; gamma=kappa*u_star*wall_distance; incoming reference strain; same diffusion and fields")},
        {"semantics",J::Quote("refreshed same-stage tensors; physical gradients on identified native elements; no averaging, lifting or acceptance waiver")},
        {"governing_solves","0"},{"physical_commits","0"},{"selected_state_unchanged","true"}});
}

std::string RANSWallEvolutionPatch::InterfaceProvenanceChecks() const
{
    namespace J=RANSWallSolveJson;
    Require(AcceptedAWF() && !Mapped() && m_dimension==2,"fitted-interface fixture requires accepted flat2D");
    const double H=m_parameters.geometry.height,K0=m_parameters.model.Q()*m_robustUstar*m_robustUstar;
    std::vector<RANSWallEvolutionFunction> functions{
        [H](const RANSWallCoordinate &p,double){return (p[0]<=-3?.01:.02)*(p[0]+3)+(p[2]<=H?.1*p[2]:.1*H+.3*(p[2]-H));},
        [](const RANSWallCoordinate&,double){return 0.;},
        [H,K0](const RANSWallCoordinate &p,double){return K0+(p[2]<=H?.01:.02)*(p[2]-H);},
        [](const RANSWallCoordinate&,double){return .1;},
        [](const RANSWallCoordinate&,double){return 0.;},
        [](const RANSWallCoordinate&,double){return 0.;}};
    const auto state=Project(functions,0.);const auto hash=RANSWallR3::Hash(state);
    RANSWallEvolutionData data(2);data.fixedBoundary=FixedReferenceBoundary(state);
    const auto domain=ScalarDomain(state,0.,data,true);
    Require(domain.admissible,"interface fixture scalar domain: "+domain.failure);
    int statics=0,dynamics=0,outer=0;double error=0;std::vector<std::string> rows;
    InitialiseScalarLocations();
    for(const auto &p:m_scalarLocations)if(p.inventory=="physical_sampled_K")
    {ScalarObservation q;q.coordinate=p.row.coordinate;q.block=p.block;q.geometry=p.row.geometry;q.side=p.side;
     const auto v=ResolveObservation(q);Require(q.side==-1 && v.inner,"static endpoint sampler is not H-minus");++statics;}
    for(const auto &p:domain.observations)
    {
        if(p.inventory=="physical_sampled_K")
        {const auto v=ResolveObservation(p);Require(p.side==-1 && v.inner,"dynamic endpoint sampler is not H-minus");++dynamics;}
        if(p.inventory=="constitutive_matching_K")
        {const auto v=ResolveObservation(p);Require(p.side==1 && !v.inner,"explicit H-plus matching observation lost");++outer;}
    }
    int neighbors=0,rejected=0;
    for(int e=0;e<fields[0]->Elements().size();++e)
    {
        const auto &el=fields[0]->Elements()[e];
        if(-3<el.lower[0]-1e-12 || -3>el.upper[0]+1e-12 || H<el.lower[2]-1e-12 || H>el.upper[2]+1e-12)continue;
        const bool inner=el.upper[2]<=H+1e-12;
        ScalarObservation q;q.coordinate={-3,0,H};q.block=K();q.geometry=el.geometry;q.side=inner?-1:1;
        const auto provenance=ResolveObservation(q);const auto &v=*provenance.velocity;
        Require(v.element==e && v.geometry==el.geometry && provenance.inner==inner,"native one-sided identity changed");
        const auto gradient=v.Gradient(state,offset[0]);
        const double gx=el.upper[0]<=-3+1e-12?.01:.02,gz=inner?.1:.3;
        const auto kr=GeometryPoint(K(),el.geometry,q.coordinate);
        const double value=kr.Value(state,offset[K()]);
        error=std::max({error,std::abs(gradient[0]-gx),std::abs(gradient[2]-gz),std::abs(value-K0),std::abs(v.Value(state)-.1*H)});
        rows.push_back(J::Object({{"geometry",J::Number(v.geometry)},{"element",J::Number(v.element)},
            {"side",J::Number(provenance.radialSide)},{"inner",J::Boolean(inner)},
            {"gradient",J::Numbers(gradient)},{"k",J::Number(value)}}));
        q.side=-q.side;try{ResolveObservation(q);}catch(const std::exception&){++rejected;}
        ++neighbors;
    }
    Require(neighbors==4 && rejected==4 && statics>0 && dynamics>=statics && outer>0 && error<2e-10,
        "fitted-interface polynomial/provenance fixture failed");
    Require(hash==RANSWallR3::Hash(state) && InitializationWork()==Vector({0.,0.}),"interface fixture changed fields or called Lift");
    return J::Object({{"audit_revision",J::Quote(AuditRevision())},{"passed","true"},
        {"static_Hminus_samples",J::Number(statics)},{"dynamic_Hminus_samples",J::Number(dynamics)},
        {"Hplus_matching_samples",J::Number(outer)},{"tangential_radial_neighbors",J::List(rows)},
        {"incompatible_side_rejections",J::Number(rejected)},{"polynomial_error",J::Number(error)},
        {"state_sha256",J::Quote(hash)},{"governing_solves","0"},{"physical_commits","0"},{"initialization_lift_calls","0"}});
}

std::string RANSWallEvolutionPatch::CompareOuterOperators(const Vector &s,double time,const RANSWallEvolutionData &data) const
{
    namespace J=RANSWallSolveJson;const auto before=RANSWallR3::Hash(s);
    Require(ScalarDomain(s,time,data).admissible,"operator counterfactual outside scalar/budget domain");
    const auto selected=ResidualResolutionComponents(s,s,0,time,data,false,false);
    const auto original=ResidualResolutionComponents(s,s,0,time,data,false,true);
    const auto &a=selected.terms.at("total_constrained"),&b=original.terms.at("total_constrained");
    Vector delta(Size());for(int i=0;i<Size();++i)delta[i]=a[i]-b[i];
    const auto native=Residual(s,s,0,time,data,false,true);double nativeDifference=0,otherDifference=0;
    for(int i=0;i<Size();++i)nativeDifference=std::max(nativeDifference,std::abs(a[i]-native[i]));
    Require(nativeDifference<1e-9,"selected assembled operator differs from full native residual");
    for(const auto &name:RANSWallResolutionSpatialTerms())
        if(name!="outer_momentum_stress" && name!="outer_k_production" && name!="epsilon_production")
            for(int i=0;i<Size();++i)otherDifference=std::max(otherDifference,std::abs(selected.terms.at(name)[i]-original.terms.at(name)[i]));
    Require(otherDifference==0.,"outer counterfactual changed an unrelated physical term");
    if(m_variantName=="Original")Require(RANSWallAlgebra::Maximum(delta)==0.,"Original same-state operator counterfactual is not an exact no-op");
    const double us=m_robustUstar,H=m_parameters.geometry.height,K0=m_parameters.model.Q()*us*us;
    std::vector<std::string> outputs;const auto na=Norms(a,selected.constraints,us),nb=Norms(b,selected.constraints,us),nd=Norms(delta,selected.constraints,us);
    for(int block=0;block<InnerPressure();++block)
    {
        Vector rows,la,lb,ld;for(int i=0;i<size[block];++i)if(!selected.constraints.count(offset[block]+i))
        {rows.push_back(offset[block]+i);la.push_back(a[offset[block]+i]);lb.push_back(b[offset[block]+i]);ld.push_back(delta[offset[block]+i]);}
        const double scale=block<m_dimension?us:(block==K()?K0:us*us*us/(m_parameters.model.Kappa()*H));
        const double factor=scale*us/H*std::sqrt(m_volumes[block]);
        const double common=(block<m_dimension?us:(block==K()?K0:us*us*us/(m_parameters.model.Kappa()*2)))*us/2*std::sqrt(block==Epsilon()?12.:24.);
        outputs.push_back(J::Object({{"block",J::Number(block)},{"free_rows",J::Numbers(rows)},{"selected_load",J::Numbers(la)},
            {"same_state_Original_load",J::Numbers(lb)},{"difference_load",J::Numbers(ld)},
            {"mass_dual_dimensional",J::Numbers(Vector{na[block]*factor,nb[block]*factor,nd[block]*factor})},
            {"common_normalization",J::Number(common)},{"common_normalized",J::Numbers(Vector{na[block]*factor/common,nb[block]*factor/common,nd[block]*factor/common})}}));
    }
    Vector constraintsA, constraintsB;for(int i=offset[InnerPressure()];i<Size();++i){constraintsA.push_back(a[i]);constraintsB.push_back(b[i]);}
    Require(before==RANSWallR3::Hash(s),"operator diagnostic changed retained fields");
    return J::Object({{"state_sha256",J::Quote(before)},{"selected_outer",J::Quote(m_variantName)},
        {"counterfactual_outer",J::Quote("Original")},{"alpha","0"},{"blocks",J::List(outputs)},
        {"selected_native_residual_maximum_difference",J::Number(nativeDifference)},
        {"unaffected_physical_term_maximum_difference",J::Number(otherDifference)},
        {"selected_continuity_gauge_loads",J::Numbers(constraintsA)},{"counterfactual_continuity_gauge_loads",J::Numbers(constraintsB)},
        {"mass_factor_applications",J::Number(3*Gauge())},{"governing_solves","0"},{"native_history_advanced","false"},
        {"legacy_name_semantics",J::Quote("refreshed_original_block_norms means selected uncancelled operator, not Original-CLS counterfactual")}});
}

std::string RANSWallEvolutionPatch::EvaluationReuseChecks(const Vector &state,const Vector &history,
    double time,const RANSWallEvolutionData &data)
{
    namespace J=RANSWallSolveJson;
    Require(!m_evaluation,"reuse verification needs an independent scope");
    const bool saved=m_evaluationReuse;const auto stateHash=RANSWallR3::Hash(state),historyHash=RANSWallR3::Hash(history);
    struct Restore {bool &value;bool old;~Restore(){value=old;}} restore{m_evaluationReuse,saved};
    const auto drive=[&](bool enabled)
    {
        m_evaluationReuse=enabled;RANSWallEvolutionData input=data;
        EvaluationScope scope(*this,input,time);std::vector<std::string> values;
        const auto domain=[&](const Vector &s,double t,bool retain,bool only,const std::vector<RANSWallCoordinate> &anchors=std::vector<RANSWallCoordinate>{})
        {values.push_back(ScalarDomainJSON(ScalarDomain(s,t,input,retain,only,anchors),retain));};
        domain(state,time,true,true);domain(state,time,true,false);domain(state,time,true,false);
        domain(state,time,false,false);
        auto changed=state;auto zero=std::find(changed.begin(),changed.end(),0.);
        Require(zero!=changed.end(),"reuse fixture needs an exact zero coefficient");*zero=std::signbit(*zero)?0.:-0.;
        domain(changed,time,true,false);domain(changed,time,true,false);
        auto flow=state,kstate=state,eps=state;
        for(const auto &[block,increment]:std::array<std::pair<int,double>,3>{{{0,1e-7},{K(),1e-7},{Epsilon(),1e-8}}})
        {
            auto &test=block==0?flow:block==K()?kstate:eps;bool found=false;
            for(int j=0;j<size[block];++j)if(!input.fixedBoundary.count(offset[block]+j))
            {test[offset[block]+j]+=increment;found=true;break;}
            Require(found,"reuse verification needs a free field coefficient");domain(test,time,true,false);
        }
        domain(flow,time,true,false);

        domain(state,std::nextafter(time,1e300),true,false);domain(state,time,true,false);
        ++input.generation;domain(state,time,true,false);
        auto anchors=std::vector<RANSWallCoordinate>{{(lower[0]+upper[0])*.5+.0123456789,0,0}};
        domain(state,time,true,true,anchors);domain(state,time,true,true,anchors);
        domain(state,time,true,false);
        {CoverageScope coverage(m_assertKCoverage,true);domain(state,time,true,false);domain(state,time,true,false);}
        for(const auto &[axis,position,side]:std::array<std::tuple<int,double,int>,6>{{{2,0.,-1},{2,-0.,-1},
            {2,m_parameters.geometry.height,-1},{2,m_parameters.geometry.height,1},{0,-3.,-1},{0,-3.,1}}})
        {
            const auto faces=Surface(0,axis,position,side,8);std::vector<std::string> rows;
            for(const auto &face:faces)rows.push_back(J::Object({{"coordinate",J::Numbers(face.test.coordinate)},
                {"geometry",J::Number(face.test.geometry)},{"element",J::Number(face.test.element)},
                {"phi",J::Numbers(face.test.phi)},{"derivative",J::Rows(face.test.derivative)}}));
            values.push_back(J::List(rows));
        }
        auto foreign=input;foreign.Touch();
        values.push_back(ScalarDomainJSON(ScalarDomain(state,time,foreign,true),true));
        domain(state,time,true,false);
        if(!input.fixedBoundary.empty())
        {input.fixedBoundary.begin()->second+=1e-7;domain(state,time,true,false);input=data;input.Touch();}
        for(int repeat=0;repeat<2;++repeat)values.push_back(J::Numbers(Residual(state,history,1600.,time,input)));
        auto h=history;h[0]+=1e-5;
        values.push_back(J::Numbers(Residual(state,h,1600.,time,input)));
        values.push_back(J::Numbers(Residual(state,history,800.,time,input,false,false,0)));
        input.innerKSource=[](const RANSWallCoordinate&,double){return .001;};input.Touch();
        values.push_back(J::Numbers(Residual(state,history,1600.,time,input)));
        std::vector<std::pair<std::string,std::string>> audit;
        for(const auto &[key,value]:CompleteRestrictedAudit(state,history,1600.,time,input))audit.push_back({key,J::Numbers(value)});
        values.push_back(J::Object(audit));
        values.push_back(CheckConstitutiveCoverage(state,time,input,{lower[0]+.123456789,0,m_parameters.geometry.height*.314159}));
        return values;
    };
    const auto off=drive(false),on=drive(true);
    Require(off.size()==on.size(),"evaluation verification query count differs");
    for(int i=0;i<off.size();++i)Require(off[i]==on[i],"exact reuse changed query "+std::to_string(i));
    Require(stateHash==RANSWallR3::Hash(state) && historyHash==RANSWallR3::Hash(history),"evaluation verification modified inputs");
    Require(m_evaluationWork.domainHits>0 && m_evaluationWork.prefixHits>0 && m_evaluationWork.inputInvalidations>0 &&
        m_evaluationWork.residualHits>0,"verification did not exercise cache hits and invalidation");
    return J::Object({{"passed","true"},{"audit_revision",J::Quote(AuditRevision())},{"compared_queries",J::Number(off.size())},
        {"state_sha256",J::Quote(stateHash)},{"history_sha256",J::Quote(historyHash)},
        {"checks",J::Quote("exact state bytes including signed zero; time; data generation; external data; flags; ordered additional anchors; foreign data; history; alpha; residual options; coverage tripwire")},
        {"reuse",EvaluationReuseJSON()},{"governing_solves","0"},{"physical_commits","0"}});
}

}
