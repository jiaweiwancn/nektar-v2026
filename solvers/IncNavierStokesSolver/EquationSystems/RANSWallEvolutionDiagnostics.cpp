///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallTiming.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionResolution.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
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
// Incompressible rank-one momentum action for propagation in the physical x-z
// plane. Reflection symmetry decouples the transverse polarization in this 2D
// pilot. A coarse bracket followed by bounded minimization is a numerical
// diagnostic on this plane, not a global 3D ellipticity certificate.
struct WJPrincipalResult { std::array<double,2> value,angle; };
WJPrincipalResult WJPrincipal(double k,double epsilon,const RANSTensor &g,const WJParameters &p)
{
    const auto law=EvaluateWJClosure(k,epsilon,g,p,WJGradientPolicy::DevStrainRawWorkV1);
    std::array<RANSTensor,3> h{};h[0][0][0]=1;h[0][2][2]=-1;h[1][0][2]=1;h[2][2][0]=1;
    std::array<RANSTensor,3> actions{};
    for(int j=0;j<3;++j)actions[j]=ApplyWJClosureDirection(k,epsilon,g,0,0,h[j],p,WJGradientPolicy::DevStrainRawWorkV1).reynoldsStress;
    double matrix[3][3]{};
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)for(int a=0;a<3;++a)for(int b=0;b<3;++b)
        matrix[i][j]-=h[i][a][b]*actions[j][a][b]/law.eddyViscosity;
    auto form=[&](double angle){const double c=std::cos(angle),s=std::sin(angle);const double v[3]={-c*s,-s*s,c*c};double result=0;
        for(int i=0;i<3;++i)for(int j=0;j<3;++j){result+=v[i]*matrix[i][j]*v[j];}
        return result;};
    const double pi=std::acos(-1.),spacing=pi/32;double best=1e300,bestAngle=0;
    for(int j=0;j<32;++j)
    {
        const double center=j*spacing;
        if(form(center)>form(center-spacing) || form(center)>form(center+spacing))continue;
        double lo=center-spacing,hi=center+spacing;
        for(int n=0;n<48;++n){const double x=lo+(hi-lo)*.3819660112501051,y=hi-(hi-lo)*.3819660112501051;
            if(form(x)<form(y))hi=y;else lo=x;}
        const double angle=(lo+hi)/2,value=form(angle);
        if(value<best){best=value;bestAngle=angle;}
    }
    RANSTensor hx{},hz{};hx[1][0]=1;hz[1][2]=1;
    const auto dx=ApplyWJClosureDirection(k,epsilon,g,0,0,hx,p,WJGradientPolicy::DevStrainRawWorkV1).reynoldsStress;
    const auto dz=ApplyWJClosureDirection(k,epsilon,g,0,0,hz,p,WJGradientPolicy::DevStrainRawWorkV1).reynoldsStress;
    const double a=-dx[1][0]/law.eddyViscosity,b=-(dx[1][2]+dz[1][0])/(2*law.eddyViscosity),c=-dz[1][2]/law.eddyViscosity;
    return {{best,(a+c-std::hypot(a-c,2*b))/2},
        {bestAngle,.5*std::atan2(2*b,a-c)+pi/2}};
}

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
                if(!WJ())
                {
                    const double h=1e-5*std::max(1.,std::abs(G[0][2]));auto plus=G,minus=G;plus[1][0]+=h;minus[1][0]-=h;
                    const auto cp=OuterClosure(k,epsilon,plus,v.coordinate),cm=OuterClosure(k,epsilon,minus,v.coordinate);
                    r[4]=std::min(r[4],-(cp.reynoldsStress[1][0]-cm.reynoldsStress[1][0])/(2*h));
                } // WJ principal diagnostics are supplied by its shared action below.
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
        double flux=0,un=0;for(int a=0;a<3;++a){flux+=((m_parameters.geometry.molecularNu+nt/SigmaK())*gki[a]-c.kDiffusion*gko[a])*n[a];}
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
namespace
{
namespace CJ=RANSWallSolveJson;
using CV=RANSWallAlgebra::Vector;
struct Spectrum { CV eigenvalues; RANSTensor vectors{}; };
Spectrum CovarianceSpectrum(RANSTensor a)
{
    RANSTensor v{};for(int i=0;i<3;++i)v[i][i]=1;
    for(int step=0;step<64;++step)
    {
        int p=0,q=1;for(int i=0;i<3;++i)for(int j=i+1;j<3;++j)if(std::abs(a[i][j])>std::abs(a[p][q])){p=i;q=j;}
        if(std::abs(a[p][q])<1e-16*std::max({1.,std::abs(a[0][0]),std::abs(a[1][1]),std::abs(a[2][2])}))break;
        const double theta=.5*std::atan2(2*a[p][q],a[q][q]-a[p][p]),c=std::cos(theta),s=std::sin(theta);
        const double pp=c*c*a[p][p]-2*s*c*a[p][q]+s*s*a[q][q],qq=s*s*a[p][p]+2*s*c*a[p][q]+c*c*a[q][q];
        for(int k=0;k<3;++k)
        {
            if(k!=p && k!=q){const double x=a[k][p],y=a[k][q];a[k][p]=a[p][k]=c*x-s*y;a[k][q]=a[q][k]=s*x+c*y;}
            const double x=v[k][p],y=v[k][q];v[k][p]=c*x-s*y;v[k][q]=s*x+c*y;
        }
        a[p][p]=pp;a[q][q]=qq;a[p][q]=a[q][p]=0;
    }
    std::array<int,3> order{0,1,2};std::sort(order.begin(),order.end(),[&](int i,int j){return a[i][i]<a[j][j];});
    Spectrum out;for(int j:order)out.eigenvalues.push_back(a[j][j]);
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)out.vectors[i][j]=v[i][order[j]];
    return out;
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
    namespace J=RANSWallSolveJson;Require(!WJ() && AcceptedAWF() && !Mapped() && m_dimension==2,"legacy constitutive forensic scope is CLS accepted AWF flat2D; use WJ native checks for WJ");
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
            c.kSource=c.production-epsilon;c.kDiffusion=m_parameters.geometry.molecularNu+nt/SigmaK();
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

namespace Nektar
{
RANSWallEvolutionPatch::Diagnostic RANSWallEvolutionPatch::WJCompleteAudit(const Vector &state,const Vector &history,
    double alpha,double time,const RANSWallEvolutionData &data,bool rangesOnly) const
{
    Require(WJ(),"WJ audit requires WJ model ownership");
    const auto domain=ScalarDomain(state,time,data,true);
    Diagnostic audit;
    audit["wj_scalar_domain"]={domain.admissible?1.:0.,domain.minimumInnerK,domain.minimumOuterK,
        domain.minimumOuterEpsilon,domain.minimumSampledK,domain.minimumTraceEpsilon,
        domain.minimumBudgetRatio,domain.maximumBudgetRelativeDefect};
    audit["wj_hard_admissible"]={0};
    audit["wj_condition_scalar_domain"]={domain.admissible?1.:0.,domain.minimumInnerK,domain.minimumOuterK,domain.minimumOuterEpsilon,domain.minimumTraceEpsilon};
    if(!domain.admissible){audit["admissible"]={0};return audit;}
    const auto scalar=audit["wj_scalar_domain"];
    if(!rangesOnly)audit=Audit(state,history,alpha,time,data,m_robustUstar,true);
    audit.erase("robust_source");audit.erase("outer_constitutive_range");
    audit["wj_scalar_domain"]=scalar;
    audit["wj_condition_scalar_domain"]={1.,scalar[1],scalar[2],scalar[3],scalar[5]};
    const double K0=m_parameters.model.Q()*m_robustUstar*m_robustUstar;
    Vector covariance{1e300,1e300},root{1e300,0,1e300},viscosity{1e300,0,1e300,0},work{0,0},principal{1e300,1e300};
    Vector invariants{1e300,-1e300,1e300,-1e300,1e300,-1e300},production{1e300,-1e300,1e300,-1e300};
    Vector innerWitness,outerWitness;int innerCount=0,outerCount=0,innerH=0,outerH=0;
    std::set<std::tuple<int,RANSWallCoordinate,double>> seen;
    for(const auto &p:domain.observations)
    {
        if(p.block!=K() && !p.dependentTrace)continue;
        const auto provenance=ResolveObservation(p);const auto &row=*provenance.velocity;
        const double k=GeometryPoint(K(),row.geometry,p.coordinate).Value(state,offset[K()]);
        RANSTensor g{};for(int a=0;a<m_dimension;++a)g[PhysicalAxis(a)]=row.Gradient(state,offset[a]);
        const bool inner=provenance.inner;
        const double epsilon=inner?0:p.dependentTrace?p.value:GeometryPoint(Epsilon(),row.geometry,p.coordinate).Value(state,offset[Epsilon()]);
        if(!seen.emplace(row.geometry,p.coordinate,epsilon).second)continue;
        const bool atH=std::abs(p.coordinate[2]-m_parameters.geometry.height)<1e-12;
        if(inner)
        {
            ++innerCount;if(atH)++innerH;
            const auto profile=StateProfile(state,p.coordinate,time,data);
            const auto law=InnerStress(k,profile->Original().EddyViscosity(Distance(p.coordinate)),g);
            const double eigen=RANSWallSolveJson::Eigenvalues(law.reynoldsStress)[0];
            if(eigen<covariance[0]){covariance[0]=eigen;innerWitness={double(row.geometry),double(row.element),double(provenance.radialSide),p.coordinate[0],p.coordinate[2],k,eigen};}
            continue;
        }
        ++outerCount;if(atH)++outerH;
        const auto law=OuterClosure(k,epsilon,g,p.coordinate);const double eigen=RANSWallSolveJson::Eigenvalues(law.reynoldsStress)[0];
        if(eigen<covariance[1]){covariance[1]=eigen;outerWitness={double(row.geometry),double(row.element),double(provenance.radialSide),p.coordinate[0],p.coordinate[2],k,epsilon,eigen};}
        root[0]=std::min(root[0],law.wj.n);root[2]=std::min(root[2],law.wj.rootDerivative);const double n=law.wj.n;
        const double iv[3]={law.wj.s,law.wj.w,law.wj.i4};
        for(int i=0;i<3;++i){invariants[2*i]=std::min(invariants[2*i],iv[i]);invariants[2*i+1]=std::max(invariants[2*i+1],iv[i]);}
        production[0]=std::min(production[0],law.production);production[1]=std::max(production[1],law.production);
        production[2]=std::min(production[2],law.productionDev);production[3]=std::max(production[3],law.productionDev);
        root[1]=std::max(root[1],std::abs(law.wj.rootResidual)/(1+n*n*n+2.7*law.wj.s*n+std::abs(2*law.wj.w*n)));
        viscosity[0]=std::min(viscosity[0],law.eddyViscosity);viscosity[1]=std::max(viscosity[1],law.eddyViscosity);
        viscosity[2]=std::min(viscosity[2],law.scalarViscosity);viscosity[3]=std::max(viscosity[3],law.scalarViscosity);
        double actual=0;for(int a=0;a<3;++a)for(int b=0;b<3;++b)actual-=law.reynoldsStress[a][b]*g[a][b];
        work[0]=std::max(work[0],std::max(std::abs(actual-law.production),std::abs(law.productionDev-law.production-2*k*law.divergence/3)));
        work[1]=std::max(work[1],std::abs(law.divergence));
        const auto values=WJPrincipal(k,epsilon,g,m_parameters.wj);
        for(int i=0;i<2;++i)if(values.value[i]<principal[i])
        {
            principal[i]=values.value[i];const double angle=values.angle[i],c=std::cos(angle),s=std::sin(angle);
            const RANSVector xi{c,0,s},a=i==0?RANSVector{-s,0,c}:RANSVector{0,1,0};RANSTensor h{};
            for(int a0=0;a0<3;++a0)for(int b=0;b<3;++b)h[a0][b]=a[a0]*xi[b];
            const auto action=ApplyWJClosureDirection(k,epsilon,g,0,0,h,m_parameters.wj,WJGradientPolicy::DevStrainRawWorkV1);
            double contraction=0,dot=0;for(int a0=0;a0<3;++a0){dot+=a[a0]*xi[a0];
                for(int b=0;b<3;++b)contraction-=h[a0][b]*action.reynoldsStress[a0][b];}
            const auto eigenvalues=RANSWallSolveJson::Eigenvalues(law.reynoldsStress);
            Vector witness{double(row.geometry),double(row.element),1.,double(provenance.radialSide),
                p.coordinate[0],p.coordinate[1],p.coordinate[2],k,epsilon,SampleK(state,p.coordinate),law.wj.n,
                eigenvalues[0],eigenvalues[1],eigenvalues[2],law.eddyViscosity,law.scalarViscosity,
                law.production,law.productionDev,law.divergence,law.wj.s,law.wj.w,law.wj.i4,law.wj.rootDerivative,
                law.wj.rootResidual,angle,principal[i],dot,contraction,contraction/law.eddyViscosity};
            witness.insert(witness.end(),xi.begin(),xi.end());witness.insert(witness.end(),a.begin(),a.end());
            for(const auto &tensor:{g,action.reynoldsStress,h})for(const auto &r:tensor)witness.insert(witness.end(),r.begin(),r.end());
            audit[i==0?"wj_principal_in_plane_witness":"wj_principal_transverse_witness"]=std::move(witness);
        }
    }
    audit["wj_normalized_invariant_s_w_I4_ranges"]=invariants;audit["wj_Praw_Pdev_ranges"]=production;
    audit["wj_covariance"]=covariance;audit["wj_root_minimum_scaled_residual"]=root;
    audit["wj_nu_eff_and_nu_d_ranges"]=viscosity;audit["wj_work_defect_raw_divergence"]=work;
    audit["wj_rank_one_xz_principal_minima_over_nu_eff"]=principal;
    audit["wj_native_element_side_counts"]={double(innerCount),double(outerCount),double(innerH),double(outerH)};
    audit["wj_inner_minimum_geometry_element_side_x_z_k_covariance"]=innerWitness;
    audit["wj_outer_minimum_geometry_element_side_x_z_k_epsilon_covariance"]=outerWitness;
    audit["wj_cls_max_beta04_R1_R2_cap_applicable"]={0,0,0,0,0};
    auto condition=[&](const std::string &name,bool passed,Vector margins)
    {margins.insert(margins.begin(),passed?1.:0.);audit["wj_condition_"+name]=std::move(margins);};
    condition("inner_covariance",covariance[0]>=-1e-10*K0,{covariance[0]+1e-10*K0});
    condition("outer_covariance",covariance[1]>=-1e-10*K0,{covariance[1]+1e-10*K0});
    condition("root_branch",root[0]>=1.8 && root[2]>0,{root[0]-1.8,root[2]});
    condition("root_residual",root[1]<1e-12,{1e-12-root[1]});
    condition("scalar_viscosities",viscosity[0]>0 && viscosity[2]>0,{viscosity[0],viscosity[2]});
    condition("principal_in_plane",principal[0]>0,{principal[0]});
    condition("principal_transverse",principal[1]>0,{principal[1]});
    condition("raw_work",work[0]<1e-11,{1e-11-work[0]});
    condition("profile_gate",scalar[6]>=1,{scalar[6]-1});
    condition("profile_budget",scalar[7]<1e-10,{1e-10-scalar[7]});
    bool valid=true;for(const auto &[name,value]:audit)if(name.find("wj_condition_")==0)valid=valid && value[0]==1;
    audit["wj_hard_admissible"]={valid?1.:0.};
    RANSTensor adverse{};adverse[0][2]=3.5;
    const auto bad=EvaluateWJClosure(1,1,adverse,m_parameters.wj);
    RANSTensor h{};const RANSVector xi{1/std::sqrt(50.),0,-7/std::sqrt(50.)},a{7/std::sqrt(50.),0,1/std::sqrt(50.)};
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)h[i][j]=a[i]*xi[j];
    const auto action=ApplyWJClosureDirection(1,1,adverse,0,0,h,m_parameters.wj);
    double form=0;for(int i=0;i<3;++i)for(int j=0;j<3;++j)form-=h[i][j]*action.reynoldsStress[i][j]/bad.eddyViscosity;
    RANSTensor commuting{{{2,1.5*std::sqrt(17.),0},{-1.5*std::sqrt(17.),2,0},{0,0,-4}}};
    const auto psd=EvaluateWJClosure(1,1,commuting,m_parameters.wj);
    audit["wj_separate_analytic_adverse_witnesses"]={bad.wj.n,form,RANSWallSolveJson::Eigenvalues(bad.reynoldsStress)[0],bad.production,
        psd.wj.n,psd.reynoldsStress[0][0],psd.reynoldsStress[1][1],psd.reynoldsStress[2][2],psd.production};
    return audit;
}
}

namespace Nektar
{
std::string RANSWallEvolutionPatch::WJAuditPacket(const Diagnostic &audit,const Vector &state,const Vector &history,
    double alpha,double time,const std::string &phase) const
{
    namespace J=RANSWallSolveJson;std::vector<std::pair<std::string,std::string>> conditions,witnesses;
    std::vector<std::string> failures;
    for(const auto &[name,v]:audit)if(name.find("wj_condition_")==0)
    {
        const std::string category=name.substr(13);const bool pass=v[0]==1;
        conditions.push_back({category,J::Object({{"passed",J::Boolean(pass)},
            {"margins",J::Numbers(Vector(v.begin()+1,v.end()))}})});
        if(!pass)failures.push_back(J::Quote(category));
    }
    for(const std::string name:{"scalar_domain","inner_covariance","outer_covariance","root_branch","root_residual",
        "scalar_viscosities","principal_in_plane","principal_transverse","raw_work","profile_gate","profile_budget"})
        if(!audit.count("wj_condition_"+name))conditions.push_back({name,J::Object({{"passed","null"},{"margins","[]"},{"evaluation",J::Quote("not_evaluated_after_scalar_domain_exit")}})});
    for(const std::string polarization:{"in_plane","transverse"})
    {
        const auto found=audit.find("wj_principal_"+polarization+"_witness");if(found==audit.end())continue;
        const auto &v=found->second;Require(v.size()==62,"WJ principal witness layout mismatch");
        auto range=[&](int a,int b){return J::Numbers(Vector(v.begin()+a,v.begin()+b));};
        witnesses.push_back({polarization,J::Object({{"geometry",J::Number(v[0])},{"element",J::Number(v[1])},
            {"region",J::Quote("outer")},{"side",J::Number(v[3])},{"position",range(4,7)},
            {"k",J::Number(v[7])},{"epsilon",J::Number(v[8])},{"sampled_K",J::Number(v[9])},
            {"N",J::Number(v[10])},{"covariance_eigenvalues",range(11,14)},{"nu_eff",J::Number(v[14])},
            {"nu_d",J::Number(v[15])},{"Praw",J::Number(v[16])},{"Pdev",J::Number(v[17])},
            {"raw_divergence",J::Number(v[18])},{"s_w_I4",range(19,22)},{"F_N",J::Number(v[22])},
            {"root_residual",J::Number(v[23])},{"angle_radians",J::Number(v[24])},
            {"minimum_over_nu_eff",J::Number(v[25])},{"a_dot_xi",J::Number(v[26])},
            {"minus_H_contract_dR",J::Number(v[27])},{"replayed_over_nu_eff",J::Number(v[28])},
            {"xi",range(29,32)},{"polarization",range(32,35)},{"raw_G_row_major",range(35,44)},
            {"dR_row_major",range(44,53)},{"H_row_major",range(53,62)},
            {"normalization",J::Quote("unit xi, unit a perpendicular to xi; -H:dR/nu_eff; fixed positive k and epsilon")},
            {"search",J::Quote(polarization=="in_plane"?"32 periodic angular brackets; 48 golden iterations per local minimum":"analytic symmetric 2x2 eigenvalue in x-z propagation plane")},
            {"scope",J::Quote("finite NativeElementSideV2 spatial inventory; x-z propagation only")}})});
    }
    return J::Object({{"schema",J::Quote("WJFiniteAuditV2")},{"state_sha256",J::Quote(RANSWallR3::Hash(state))},
        {"stage_context_sha256",J::Quote(RANSWallR3::Hash(RANSWallR3::Hash(state)+RANSWallR3::Hash(history)+J::Number(alpha)+"/"+J::Number(time)+"/"+phase))},
        {"stage_history_sha256",J::Quote(RANSWallR3::Hash(history))},{"callback_history",alpha>0?J::Numbers(history):"null"},{"phase",J::Quote(phase)},
        {"time",J::Number(time)},{"alpha_mass",J::Number(alpha)},{"conditions",J::Object(conditions)},
        {"failure_categories",J::List(failures)},{"principal_witnesses",J::Object(witnesses)}});
}
}

namespace Nektar
{
std::string RANSWallEvolutionPatch::WJMeasurements(const Vector &state,const Vector &stageHistory,double stageAlpha,double time,
    const RANSWallEvolutionData &data,const std::vector<std::array<double,5>> &observations) const
{
    namespace J=RANSWallSolveJson;
    Require(WJ() && m_eta==0 && stageAlpha>=0,"WJ measurement requires unanchored selected model");
    const auto stateHash=RANSWallR3::Hash(state),historyHash=RANSWallR3::Hash(stageHistory),referenceHash=RANSWallR3::Hash(m_reference);
    const auto constraints=Constraints(state,time,data);
    const auto spatial=Residual(state,state,0,time,data,false);
    const auto resolution=ResidualResolutionComponents(state,state,0,time,data,false);
    std::vector<std::pair<std::string,std::string>> terms;
    for(const auto &[name,value]:resolution.terms)terms.push_back({name,J::Object({{"load",J::Numbers(value)},
        {"mass_dual_norms",J::Numbers(Norms(value,constraints,m_robustUstar))}})});
    std::vector<std::string> samples,interfaces;double divergence=0;Vector evaluationSensitivity(6);int observed=0;
    for(const auto &o:observations)
    {
        const RANSWallCoordinate p{o[0],o[1],o[2]};const int side=int(o[4]);
        Require(p[2]!=m_parameters.geometry.height || side!=0,"WJ H observation must name a side");
        const bool inner=p[2]<m_parameters.geometry.height || (p[2]==m_parameters.geometry.height && side<0);
        const auto v=Point(0,p,side),pr=Point(inner?InnerPressure():OuterPressure(),p,side);
        RANSTensor g{};for(int a=0;a<m_dimension;++a)g[PhysicalAxis(a)]=v.Gradient(state,offset[a]);
        divergence=std::max(divergence,std::abs(g[0][0]+g[2][2]));
        auto values=[&](const Vector &s){const double k=v.Value(s,offset[K()]),pi=pr.Value(s,offset[inner?InnerPressure():OuterPressure()]);
            return Vector{v.Value(s,offset[0]),v.Value(s,offset[1]),k,inner?0:Point(Epsilon(),p,side).Value(s,offset[Epsilon()]),pi,WJMechanicalPressure(pi,k,0)};};
        const auto value=values(state),initial=values(m_reference);
        const int blocks[5]={0,1,K(),Epsilon(),inner?InnerPressure():OuterPressure()};Vector alternate(6);
        if(observed++%std::max(1,int(observations.size()/64))==0)
        {for(int j=0;j<5;++j)if(j!=3 || !inner)alternate[j]=fields[blocks[j]]->NativeEvaluate(state,p,offset[blocks[j]],side)[0];
        alternate[5]=WJMechanicalPressure(alternate[4],alternate[2],0);
        for(int j=0;j<6;++j)evaluationSensitivity[j]=std::max(evaluationSensitivity[j],std::abs(value[j]-alternate[j]));}
        Vector G;for(const auto &r:g)G.insert(G.end(),r.begin(),r.end());
        samples.push_back(J::Object({{"position",J::Numbers(p)},{"weight",J::Number(o[3])},{"side",J::Number(side)},
            {"inner",J::Boolean(inner)},{"value",J::Numbers(value)},{"represented_initial",J::Numbers(initial)},
            {"raw_G",J::Numbers(G)},{"grad_k",J::Numbers(v.Gradient(state,offset[K()]))}}));
    }
    for(const auto &f:Surface(0,2,m_parameters.geometry.height,1,8))
    {
        const auto p=f.test.coordinate;const auto vi=Point(0,p,-1),vo=Point(0,p,1);
        RANSTensor gi{},go{};for(int a=0;a<m_dimension;++a){gi[PhysicalAxis(a)]=vi.Gradient(state,offset[a]);go[PhysicalAxis(a)]=vo.Gradient(state,offset[a]);}
        const double ki=vi.Value(state,offset[K()]),ko=vo.Value(state,offset[K()]);
        const double epsilon=Point(Epsilon(),p,1).Value(state,offset[Epsilon()]);
        const double pii=Point(InnerPressure(),p,-1).Value(state,offset[InnerPressure()]),pio=Point(OuterPressure(),p,1).Value(state,offset[OuterPressure()]);
        const auto profile=StateProfile(state,p,time,data);const double nt=profile->Original().EddyViscosity(p[2]);
        const auto ti=RANSWallInnerMeanStress(InnerStress(ki,nt,gi),m_parameters.geometry.molecularNu,pii);
        const auto co=OuterClosure(ko,epsilon,go,p);
        Vector innerTraction(3),outerTraction(3);for(int a=0;a<3;++a){innerTraction[a]=ti[a][2];outerTraction[a]=co.reynoldsStress[a][2]-(a==2?2*ko/3-pio:0);}
        interfaces.push_back(J::Object({{"position",J::Numbers(p)},{"weight",J::Number(f.weight)},
            {"inner_outer_k",J::Numbers(Vector{ki,ko})},{"epsilon_trace_target",J::Numbers(Vector{epsilon,profile->Dissipation(p[2])})},
            {"sampled_K",J::Number(profile->K())},{"cutoff",J::Number(profile->Reference().cutoff)},
            {"inner_outer_pi",J::Numbers(Vector{pii,pio})},
            {"inner_outer_mechanical_pressure",J::Numbers(Vector{WJMechanicalPressure(pii,ki,0),WJMechanicalPressure(pio,ko,0)})},
            {"inner_outward_traction",J::Numbers(innerTraction)},{"outer_outward_traction",J::Numbers(outerTraction)},
            {"opposing_k_diffusive_flux",J::Numbers(Vector{(m_parameters.geometry.molecularNu+nt/SigmaK())*vi.Gradient(state,offset[K()])[2],-co.kDiffusion*vo.Gradient(state,offset[K()])[2]})}}));
    }
    const auto audit=WJCompleteAudit(state,stageAlpha>0?stageHistory:state,stageAlpha,time,data,true);
    const auto stage=stageAlpha>0?Residual(state,stageHistory,stageAlpha,time,data,false):Vector{};
    Require(stateHash==RANSWallR3::Hash(state) && historyHash==RANSWallR3::Hash(stageHistory) && referenceHash==RANSWallR3::Hash(m_reference) && InitializationWork()==Vector({0.,0.}),"WJ observation mutated immutable state/history/reference or called Lift");
    return J::Object({{"schema",J::Quote("WJReadOnlyComparisonV1")},{"state_sha256",J::Quote(stateHash)},
        {"state_after_sha256",J::Quote(RANSWallR3::Hash(state))},{"stage_history_sha256",J::Quote(historyHash)},
        {"stage_history_after_sha256",J::Quote(RANSWallR3::Hash(stageHistory))},{"represented_initial_sha256",J::Quote(referenceHash)},
        {"time",J::Number(time)},{"alpha_mass_spatial","0"},{"eta","0"},{"governing_solves","0"},{"initialization_lift_calls_iterations",J::Numbers(InitializationWork())},
        {"spatial_residual",J::Numbers(spatial)},{"spatial_mass_dual_norms",J::Numbers(Norms(spatial,constraints,m_robustUstar))},
        {"completed_stage",stageAlpha>0?J::Object({{"alpha_mass",J::Number(stageAlpha)},{"time",J::Number(time)},
            {"history_sha256",J::Quote(historyHash)},{"residual",J::Numbers(stage)},{"mass_dual_norms",J::Numbers(Norms(stage,constraints,m_robustUstar))}}):"null"},
        {"regional_spatial_mass_dual_norms",J::Rows(std::vector<Vector>{Norms(Residual(state,state,0,time,data,false,true,0),constraints,m_robustUstar),Norms(Residual(state,state,0,time,data,false,true,1),constraints,m_robustUstar)})},
        {"spatial_terms",J::Object(terms)},{"constraints",J::Constraints(constraints)},
        {"residual_order",J::Quote("u,v,k,entire outer epsilon,inner continuity,outer continuity,gauge,boundary/trace; consistent mass dual physical scaling")},
        {"component_norm_warning",J::Quote("separately normed regional/source/diffusion/interface contributions are not linearly additive")},
        {"observations",J::List(samples)},{"interfaces",J::List(interfaces)},
        {"strong_divergence_sampled_Linf",J::Number(divergence)},{"arithmetic_evaluation_sensitivity",J::Numbers(evaluationSensitivity)},
        {"audit",WJAuditPacket(audit,state,stageAlpha>0?stageHistory:state,stageAlpha,time,stageAlpha>0?"accepted_checkpoint_stage":"startup_readback")}});
}
}
