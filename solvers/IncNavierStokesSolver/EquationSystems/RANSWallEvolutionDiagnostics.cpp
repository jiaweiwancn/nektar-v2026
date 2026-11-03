///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <algorithm>
#include <cmath>
#include <limits>

namespace Nektar
{
using namespace SolverUtils;
namespace
{
double MinimumEigenvalue(RANSTensor a)
{
    // Symmetric Jacobi rotations; covariance is a real symmetric tensor.
    for(int step=0;step<32;++step)
    {
        int p=0,q=1;for(int i=0;i<3;++i)for(int j=i+1;j<3;++j)
            if(std::abs(a[i][j])>std::abs(a[p][q])){p=i;q=j;}
        if(std::abs(a[p][q])<1e-15*std::max({1.,std::abs(a[0][0]),std::abs(a[1][1]),std::abs(a[2][2])}))break;
        const double angle=.5*std::atan2(2*a[p][q],a[q][q]-a[p][p]),c=std::cos(angle),s=std::sin(angle);
        const double pp=c*c*a[p][p]-2*s*c*a[p][q]+s*s*a[q][q];
        const double qq=s*s*a[p][p]+2*s*c*a[p][q]+c*c*a[q][q];
        for(int k=0;k<3;++k)if(k!=p && k!=q)
        {const double x=a[k][p],y=a[k][q];a[k][p]=a[p][k]=c*x-s*y;a[k][q]=a[q][k]=s*x+c*y;}
        a[p][p]=pp;a[q][q]=qq;a[p][q]=a[q][p]=0;
    }
    return std::min({a[0][0],a[1][1],a[2][2]});
}
}
RANSWallEvolutionPatch::Diagnostic RANSWallEvolutionPatch::Audit(const Vector &s,const Vector &history,
    double alpha,double time,const RANSWallEvolutionData &data) const
{
    Diagnostic result;const double sentinel=std::numeric_limits<double>::max();
    bool admissible=true;Refresh(s,time,&data);
    for(const std::string name:{"source_volume_inner","source_volume_outer","physical_sample_K","physical_sample_outer",
        "prescribed_source_inner","mass_volume","interface_inner","interface_outer","reconstruction_dense","reconstruction_native"})
        result[name]={0,sentinel,sentinel,sentinel,sentinel,0};
    // Every rank evaluates replicated diagnostic sets. Owned source sets are
    // reduced below; replicated sets are counted only on rank zero.
    auto sample=[&](const std::string &name,const RANSWallGlobalPoint &v,bool inner,bool constitutive,bool owned)
    {
        if(!owned && m_comm->GetRank()!=0)return;
        auto &r=result[name];if(r.empty())r={0,sentinel,sentinel,sentinel,sentinel,0};
        ++r[0];const double k=v.Value(s,offset[K()]);r[1]=std::min(r[1],k);
        const double epsilon=inner?(constitutive?DependentEpsilon(s,v.coordinate,time,data):sentinel):
            Point(Epsilon(),v.coordinate,1).Value(s,offset[Epsilon()]);
        r[2]=std::min(r[2],epsilon);
        if(!std::isfinite(k) || k<0 || (!inner && !(k>0)) || !std::isfinite(epsilon) || !(epsilon>0))++r[5];
        if(constitutive && k>0 && epsilon>0 && std::isfinite(k) && std::isfinite(epsilon))
        {
            RANSTensor G{},covariance{};for(int a=0;a<m_dimension;++a)G[PhysicalAxis(a)]=v.Gradient(s,offset[a]);
            if(inner)
            {
                const double nt=StateProfile(s,v.coordinate,time,data)->Original().EddyViscosity(Distance(v.coordinate));
                for(int a=0;a<3;++a)for(int b=0;b<3;++b)covariance[a][b]=(a==b?2*k/3:0)-nt*(G[a][b]+G[b][a]);
            }
            else
            {
                const auto c=EvaluateRANSClosure(m_parameters.model.Model(),k,epsilon,G,Distance(v.coordinate),m_parameters.model.Parameters(),0);
                covariance=c.reynoldsStress;
                const double h=1e-5*std::max(1.,std::abs(G[0][2]));auto plus=G,minus=G;plus[1][0]+=h;minus[1][0]-=h;
                const auto cp=EvaluateRANSClosure(m_parameters.model.Model(),k,epsilon,plus,Distance(v.coordinate),m_parameters.model.Parameters(),0);
                const auto cm=EvaluateRANSClosure(m_parameters.model.Model(),k,epsilon,minus,Distance(v.coordinate),m_parameters.model.Parameters(),0);
                r[4]=std::min(r[4],-(cp.reynoldsStress[1][0]-cm.reynoldsStress[1][0])/(2*h));
            }
            r[3]=std::min(r[3],MinimumEigenvalue(covariance));
        }
    };
    for(const auto &ray:m_rays)
    {
        for(const auto &q:ray.points)sample("source_volume_inner",q.v,true,true,true);
        sample("physical_sample_K",ray.sample,true,true,false);
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
    for(int e=0;e<fields[0]->Elements().size();++e)
    {
        const auto &el=fields[0]->Elements()[e];const auto exp=fields[0]->Field()->GetExp(e);
        for(int kind=0;kind<2;++kind)
        {
            const int n=kind?2*el.degree+3:exp->GetNumPoints(0);
            for(int k=0;k<(m_dimension==3?n:1);++k)for(int j=0;j<n;++j)for(int i=0;i<n;++i)
            {
                RANSWallCoordinate p{};const std::array<int,3> index{i,m_dimension==3?j:0,m_dimension==3?k:j};
                for(int a=0;a<3;++a)if(m_dimension==3 || a!=1)
                {
                    const int nativeAxis=m_dimension==2 && a==2?1:a;
                    const double x=kind?-1.+2.*index[a]/(n-1):exp->GetBasis(nativeAxis)->GetZ()[index[a]];
                    p[a]=el.lower[a]+(x+1)*(el.upper[a]-el.lower[a])/2;
                }
                sample(kind?"reconstruction_dense":"reconstruction_native",Transform(fields[0]->ElementPoint(e,p)),
                    el.upper[2]<=m_parameters.geometry.height,false,false);
            }
        }
    }
    for(auto &[name,r]:result)
    {
        (void)name;for(int i:{0,5})m_comm->AllReduce(r[i],LibUtilities::ReduceSum);
        for(int i=1;i<=4;++i)m_comm->AllReduce(r[i],LibUtilities::ReduceMin);
        if(r[5])admissible=false;
        // Finite JSON placeholder for quantities not evaluated on this set.
        for(int i=1;i<=4;++i)if(r[i]==sentinel)r[i]=0;
    }
    result["admissible"]={admissible?1.:0.};
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
        const auto c=EvaluateRANSClosure(m_parameters.model.Model(),k,eps,go,Distance(p),m_parameters.model.Parameters(),0);
        RANSVector jump{};for(int a=0;a<3;++a)for(int b=0;b<3;++b)
            jump[a]+=((m_parameters.geometry.molecularNu+nt)*(gi[a][b]+gi[b][a])-(a==b?pi:0)+c.reynoldsStress[a][b]-(a==b?2*k/3-po:0))*n[b];
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
    result["coupled_norms"]=Norms(Residual(s,history,alpha,time,data),constraints);
    Vector balance(2*(m_dimension+1));
    for(int b=0;b<=K();++b)
    {
        for(int i=0;i<size[b];++i){const int row=offset[b]+i;balance[b]+=m_constants[b][i]*raw[row];
            if(constraints.count(row))balance[m_dimension+1+b]+=m_constants[b][i]*raw[row];}
    }
    result["momentum_k_load_and_boundary_reaction"]=balance;
    Vector wall(5);for(const auto &ray:m_rays)for(const auto &q:ray.points)
    {
        RANSTensor G{};for(int a=0;a<m_dimension;++a)G[PhysicalAxis(a)]=q.v.Gradient(s,offset[a]);
        const double k=q.v.Value(s,offset[K()]),nt=ray.profile->Original().EddyViscosity(Distance(q.v.coordinate));
        double production=0,strain=0;for(int a=0;a<3;++a)for(int b=0;b<3;++b)
        {const double S=G[a][b]+G[b][a];production-=((a==b?2*k/3:0)-nt*S)*G[a][b];strain+=S*S;}
        const auto values=ray.profile->AtProduction(Distance(q.v.coordinate),production);
        wall[0]+=q.weight*ray.profile->Original().MeanOriginalDissipation();wall[1]+=q.weight*values.dissipation;
        wall[2]+=q.weight*values.redistribution;wall[3]+=q.weight*production;
        wall[4]+=q.weight*m_parameters.geometry.molecularNu*strain/2;
    }
    m_comm->AllReduce(wall,LibUtilities::ReduceSum);result["wall_sources"]=wall;
    auto energy=EnergyBudget(s,history,alpha,time,data);double reaction=0;
    for(int b=0;b<=K();++b)for(int i=0;i<size[b];++i)if(constraints.count(offset[b]+i))
        reaction+=(b==K()?m_constants[b][i]:s[offset[b]+i])*raw[offset[b]+i];
    energy.push_back(reaction);energy.push_back(AnchorEnergyWork(s));energy.push_back(energy[0]-reaction-energy[10]);
    result["energy"]=energy;result["roughness"]=RoughnessBudget(s,time,data);
    return result;
}
}
