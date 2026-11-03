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
    double alpha,double time,const RANSWallEvolutionData &data,double ustar,bool constitutiveReconstruction) const
{
    Diagnostic result;const double sentinel=std::numeric_limits<double>::max();
    Vector range{0,sentinel,-sentinel,sentinel,-sentinel,sentinel,-sentinel,0};
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
        ++r[0];const double k=v.Value(s,offset[K()]);r[1]=std::min(r[1],k);
        const double epsilon=inner?(constitutive?DependentEpsilon(s,v.coordinate,time,data):sentinel):
            Point(Epsilon(),v.coordinate,std::abs(v.coordinate[2]-m_parameters.geometry.height)<1e-12?1:0).Value(s,offset[Epsilon()]);
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
                const auto c=OuterClosure(k,epsilon,G,v.coordinate);
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
                const double h=1e-5*std::max(1.,std::abs(G[0][2]));auto plus=G,minus=G;plus[1][0]+=h;minus[1][0]-=h;
                const auto cp=OuterClosure(k,epsilon,plus,v.coordinate);
                const auto cm=OuterClosure(k,epsilon,minus,v.coordinate);
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
            const int dense=2*el.degree+(constitutiveReconstruction && el.upper[2]<=m_parameters.geometry.height?7:3);
            const std::array<int,3> counts{kind?dense:exp->GetNumPoints(0),
                m_dimension==3?(kind?dense:exp->GetNumPoints(1)):1,
                kind?dense:exp->GetNumPoints(m_dimension-1)};
            for(int k=0;k<(m_dimension==3?counts[2]:1);++k)
                for(int j=0;j<(m_dimension==3?counts[1]:counts[2]);++j)for(int i=0;i<counts[0];++i)
            {
                RANSWallCoordinate p{};const std::array<int,3> index{i,m_dimension==3?j:0,m_dimension==3?k:j};
                for(int a=0;a<3;++a)if(m_dimension==3 || a!=1)
                {
                    const int nativeAxis=m_dimension==2 && a==2?1:a;
                    const double x=kind?-1.+2.*index[a]/(counts[a]-1):exp->GetBasis(nativeAxis)->GetZ()[index[a]];
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
    if(constitutiveReconstruction)
    {
        for(int i:{0,2,4,6,7})m_comm->AllReduce(range[i],LibUtilities::ReduceMax);
        for(int i:{1,3,5})m_comm->AllReduce(range[i],LibUtilities::ReduceMin);
        result["outer_constitutive_range"]=range;
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
    m_comm->AllReduce(wall,LibUtilities::ReduceSum);result["wall_sources"]=wall;
    for(int i=0;i<6;++i)m_comm->AllReduce(repair[i],LibUtilities::ReduceSum);
    m_comm->AllReduce(repair[6],LibUtilities::ReduceMin);m_comm->AllReduce(repair[7],LibUtilities::ReduceMax);
    result["robust_source"]=repair;
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
}
