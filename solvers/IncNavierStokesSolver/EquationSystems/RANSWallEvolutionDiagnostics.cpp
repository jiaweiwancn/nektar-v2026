///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
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
}
RANSWallEvolutionPatch::Diagnostic RANSWallEvolutionPatch::Audit(const Vector &s,const Vector &history,
    double alpha,double time,const RANSWallEvolutionData &data,double ustar,bool constitutiveReconstruction) const
{
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
                const double h=1e-5*std::max(1.,std::abs(G[0][2]));auto plus=G,minus=G;plus[1][0]+=h;minus[1][0]-=h;
                const auto cp=OuterClosure(k,epsilon,plus,v.coordinate);
                const auto cm=OuterClosure(k,epsilon,minus,v.coordinate);
                r[4]=std::min(r[4],-(cp.reynoldsStress[1][0]-cm.reynoldsStress[1][0])/(2*h));
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
    namespace J=RANSWallSolveJson;Require(AcceptedAWF() && !Mapped() && m_dimension==2,"constitutive forensic scope is accepted AWF flat2D");
    const auto before=RANSWallR3::Hash(s);const auto domain=ScalarDomain(s,time,data,true);
    if(!domain.admissible)return J::Object({{"scalar_budget",ScalarDomainJSON(domain)}, {"covariance_evaluated","false"},{"state_sha256",J::Quote(before)}});
    Refresh(s,time,&data);double minimum[4]={1e300,1e300,1e300,1e300};std::string minima[4]={"null","null","null","null"};std::map<std::string,int> counts;
    auto evaluate=[&](const RANSWallGlobalPoint &v,bool inner,const std::string &family,bool retain,double epsilonOverride=std::numeric_limits<double>::quiet_NaN())
    {
        RANSTensor g{};for(int a=0;a<m_dimension;++a)g[PhysicalAxis(a)]=v.Gradient(s,offset[a]);
        const double k=v.Value(s,offset[K()]);const auto profile=StateProfile(s,v.coordinate,time,data);
        RANSWallGlobalPoint e;
        if(!inner)
        {
            e=Point(Epsilon(),v.coordinate,std::abs(v.coordinate[2]-m_parameters.geometry.height)<1e-11?1:0);
            for(int j=0;j<fields[Epsilon()]->Elements().size();++j)
                if(fields[Epsilon()]->Elements()[j].geometry==v.geometry){e=fields[Epsilon()]->ElementPoint(j,v.coordinate);break;}
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
        evaluate(RestrictedGradientRow(point),RestrictedInner(point),point.inventory,false,
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
    return J::Object({{"state_sha256",J::Quote(before)},{"scalar_budget",ScalarDomainJSON(domain)},
        {"minimum_inner",nativeInner},{"minimum_outer",nativeOuter},
        {"minimum_with_explicit_neighbors_inner",minima[0]},{"minimum_with_explicit_neighbors_outer",minima[1]},
        {"minimum_same_state_Original_outer",minima[2]},{"minimum_same_state_Regularized_outer",minima[3]},{"one_sided_witnesses",J::List(neighbors)},
        {"inventory_counts",J::Object(inventory)},{"offsets",J::Numbers(offset)},
        {"same_state_Regularized_definition",J::Quote("beta=.4; gamma=kappa*u_star*wall_distance; incoming reference strain; same diffusion and fields")},
        {"semantics",J::Quote("refreshed same-stage tensors; physical gradients on identified native elements; no averaging, lifting or acceptance waiver")},
        {"governing_solves","0"},{"physical_commits","0"},{"selected_state_unchanged","true"}});
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

}
