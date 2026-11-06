// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionResolution.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallKBoundProjection.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <cmath>
#include <set>

namespace Nektar
{
using namespace RANSWallAlgebra;
using namespace SolverUtils;
namespace J=RANSWallSolveJson;

// Diagnostic replays only. No call from a production implicit callback.
std::string RANSWallEvolutionPatch::PimpleContactDiagnostic(const Vector &s,double time,
    const RANSWallEvolutionData &data,const Vector &steps,
    const std::vector<RANSWallCoordinate> &witnesses,bool exportMass,bool isolate) const
{
    Require(PimpleSingleSweep() && ScalarDomain(s,time,data).admissible,"contact diagnostic requires an admissible PIMPLE state");
    const auto hash=RANSWallR3::Hash(s);const auto ground=WallTraceR4(s);
    const auto worst=std::min_element(ground.begin(),ground.end(),[](const auto &a,const auto &b){return a.value<b.value;});
    std::vector<RANSWallCoordinate> points{worst->row.coordinate,{1.60,0.,0.}};
    points.insert(points.end(),witnesses.begin(),witnesses.end());
    const auto constraints=Constraints(s,time,data),external=ExternalConstraints(time,data);
    auto resolution=ResidualResolutionComponents(s,s,0.,time,data,false);
    resolution.native["spatial_unmasked"]=Residual(s,s,0.,time,data,false,false);
    const auto rate=RANSWallResolutionMassRate(*this,resolution,K(),worst->row.coordinate);
    const auto rawSpatial=Residual(s,s,0.,time,data,false,true);
    double trace=0.;for(int id:m_traceIds)trace=std::max(trace,std::abs(s[offset[Epsilon()]+id]-constraints.at(offset[Epsilon()]+id)));
    const auto block=AssemblePicardBlock(K(),s,s,1600.,time,data,external);
    const auto mass=PicardKMass(block.rows);const int n=block.rows.size();
    Vector nativeRate(n),nativeLoad(n);std::vector<int> index(size[K()],-1);
    for(int i=0;i<n;++i){index[block.rows[i]-offset[K()]]=i;nativeRate[i]=rate.rate[block.rows[i]-offset[K()]];nativeLoad[i]=resolution.native.at("spatial_unmasked")[block.rows[i]];}
    double agreement=0.;for(int i=0;i<n;++i)agreement=std::max(agreement,std::abs(nativeLoad[i]-block.residual[i]));
    Require(agreement<1e-9*std::max(1.,Maximum(nativeLoad)),"contact frozen/native spatial load differs");
    auto massNorm=[&](const Vector &x){return std::sqrt(std::max(0.,RANSWallKBound::Dot(x,Multiply(mass,x))));};
    auto allNorm=[&](const Vector &x){return std::sqrt(std::max(0.,RANSWallKBound::Dot(x,Multiply(m_mass[K()],x))));};
    double nearHeight=1e300;std::set<int> nearElements;
    for(const auto &e:fields[K()]->Elements())if(e.lower[2]==0)nearHeight=std::min(nearHeight,e.upper[2]);
    for(const auto &e:fields[K()]->Elements())if(e.upper[2]<=nearHeight+1e-14)nearElements.insert(e.geometry);
    Vector integralLoad(size[K()]),nearLoad(size[K()]),fieldNorms,fieldIntegrals;
    for(int b=0;b<Gauge();++b)
    {
        double square=0,integral=0;for(const auto &q:Quadrature(b))
        {const double v=q.Value(s,offset[b]);square+=q.weight*v*v;integral+=q.weight*v;}
        fieldNorms.push_back(std::sqrt(square));fieldIntegrals.push_back(integral);
    }
    for(const auto &q:Quadrature(K()))for(int j=0;j<q.ids.size();++j)
    {integralLoad[q.ids[j]]+=q.weight*q.phi[j];if(nearElements.count(q.geometry))nearLoad[q.ids[j]]+=q.weight*q.phi[j];}
    auto integral=[&](const Vector &x,const Vector &load){return RANSWallKBound::Dot(x,load);};
    Vector kval(size[K()]);for(int i=0;i<kval.size();++i)kval[i]=s[offset[K()]+i];
    std::vector<std::string> contacts,componentIntegrals;
    for(const auto &[name,value]:rate.componentRates)
        componentIntegrals.push_back(J::Object({{"term",J::Quote(name)},{"global",J::Number(integral(value,integralLoad))},{"near",J::Number(integral(value,nearLoad))}}));
    double maximumRate=0.;Vector productionDestruction(4);
    WallMatrix reaction(n,Vector(n)),diffusion(n,Vector(n));
    auto coefficients=[&](const RANSWallGlobalPoint &q,const RANSVector &u,double w,double diff,double destruct)
    {
        (void)u;
        for(int i=0;i<q.ids.size();++i)if(index[q.ids[i]]>=0)
            for(int j=0;j<q.ids.size();++j)if(index[q.ids[j]]>=0)
            {
                reaction[index[q.ids[i]]][index[q.ids[j]]]+=w*destruct*q.phi[i]*q.phi[j];
                double dot=0.;for(int a:{0,2})dot+=q.derivative[a][i]*q.derivative[a][j];
                diffusion[index[q.ids[i]]][index[q.ids[j]]]+=w*diff*dot;
            }
    };
    for(const auto &ray:m_rays)for(const auto &q:ray.points)
    {
        RANSTensor g{};RANSVector u{};for(int a=0;a<m_dimension;++a){g[PhysicalAxis(a)]=q.v.Gradient(s,offset[a]);u[PhysicalAxis(a)]=q.v.Value(s,offset[a]);}
        const double k=q.v.Value(s,offset[K()]),nt=ray.profile->Original().EddyViscosity(Distance(q.v.coordinate));
        const auto law=InnerStress(k,nt,g);double destructionRate=0.;
        const auto src=WallSource(s,q.v.coordinate,law.production,*ray.profile,&ray.destruction,&destructionRate);
        maximumRate=std::max(maximumRate,destructionRate);
        productionDestruction[0]+=q.weight*law.production;productionDestruction[1]+=q.weight*(law.production-src.kSource);
        if(nearElements.count(q.v.geometry)){productionDestruction[2]+=q.weight*law.production;productionDestruction[3]+=q.weight*(law.production-src.kSource);}
        coefficients(q.v,u,q.weight,m_parameters.geometry.molecularNu+nt/m_parameters.model.Parameters().sigmaK,destructionRate);
    }
    for(const auto &q:m_outer)
    {
        RANSTensor g{};RANSVector u{};for(int a=0;a<m_dimension;++a){g[PhysicalAxis(a)]=q.v.Gradient(s,offset[a]);u[PhysicalAxis(a)]=q.v.Value(s,offset[a]);}
        const double k=q.v.Value(s,offset[K()]),e=q.e.Value(s,offset[Epsilon()]);const auto law=OuterClosure(k,e,g,q.v.coordinate);
        maximumRate=std::max(maximumRate,e/k);productionDestruction[0]+=q.weight*law.production;productionDestruction[1]+=q.weight*e;
        coefficients(q.v,u,q.weight,law.kDiffusion,e/k);
    }
    for(const auto &p:points)
    {
        Require(Distance(p)<=m_parameters.geometry.height,"contact witness must be in the inner region");
        const auto row=Point(K(),p);std::vector<std::pair<std::string,std::string>> terms;
        double sum=0.;for(const auto &[name,v]:rate.componentRates){const double a=row.Value(v);sum+=a;terms.push_back({name,J::Number(a)});}
        const double k=row.Value(s,offset[K()]);const auto profile=StateProfile(s,p,time,data);RANSTensor g{};
        for(int a=0;a<m_dimension;++a)g[PhysicalAxis(a)]=Point(a,p).Gradient(s,offset[a]);
        const double nt=profile->Original().EddyViscosity(Distance(p));const auto law=InnerStress(k,nt,g);
        const auto normal=WallNormalization(s,p,*profile);double c=0.;const auto source=WallSource(s,p,law.production,*profile,&normal,&c);
        contacts.push_back(J::Object({{"coordinate",J::Numbers(p)},{"geometry",J::Number(row.geometry)},{"k",J::Number(k)},
            {"mass_inverted_rate",J::Number(row.Value(rate.rate))},{"component_rates",J::Object(terms)},
            {"component_sum_error",J::Number(sum-row.Value(rate.rate))},{"point_production",J::Number(law.production)},
            {"point_net_destruction",J::Number(law.production-source.kSource)},{"point_source",J::Number(source.kSource)},
            {"point_destruction_rate",J::Number(c)},{"nu_t",J::Number(nt)}}));
    }
    std::vector<std::string> profiles;
    for(double x:{-3.,-1.,0.,1.,1.6})for(double z:{0.,.0025,.005,.01,.02,.05,.1,.25,.5,1.,2.,3.,4.,6.})if(z<=upper[2])
    {
        const RANSWallCoordinate p{x,0,z};Vector values;
        for(int b=0;b<=K();++b)values.push_back(Point(b,p).Value(s,offset[b]));
        if(z>=m_parameters.geometry.height)values.push_back(Point(Epsilon(),p,1).Value(s,offset[Epsilon()]));
        profiles.push_back(J::Object({{"coordinate",J::Numbers(p)},{"u_v_k_outer_epsilon",J::Numbers(values)}}));
    }
    auto flux=[&](double height)
    {
        Vector result(2);
        for(int a:{0,2})for(int side:{-1,1})for(const auto &face:Surface(K(),a,a==2?(side<0?0.:height):(side<0?lower[0]:upper[0]),side,16))
        {
            const auto p=face.test.coordinate;if(a==0 && p[2]>height)continue;
            const double k=face.test.Value(s,offset[K()]);RANSVector u{};RANSTensor g{};
            for(int b=0;b<m_dimension;++b){u[PhysicalAxis(b)]=Point(b,p).Value(s,offset[b]);g[PhysicalAxis(b)]=Point(b,p).Gradient(s,offset[b]);}
            double diff;
            if(p[2]<=m_parameters.geometry.height)diff=m_parameters.geometry.molecularNu+StateProfile(s,p,time,data)->Original().EddyViscosity(Distance(p))/m_parameters.model.Parameters().sigmaK;
            else diff=OuterClosure(k,Point(Epsilon(),p).Value(s,offset[Epsilon()]),g,p).kDiffusion;
            const auto grad=face.test.Gradient(s,offset[K()]);
            for(int axis:{0,2}){result[0]+=face.weight*k*u[axis]*face.area[axis];result[1]-=face.weight*diff*grad[axis]*face.area[axis];}
        }
        return result;
    };
    std::vector<std::string> scans;int solves=0;
    const std::vector<std::string> variants=isolate?std::vector<std::string>{"complete","destruction","diffusion","destruction_diffusion"}:std::vector<std::string>{"complete"};
    for(const auto &variant:variants)for(double h:steps)
    {
        Require(h>0 && std::isfinite(h),"invalid diagnostic h");WallMatrix a=mass;Vector r=nativeLoad;
        for(int i=0;i<n;++i)for(int j=0;j<n;++j)
        {
            const double l=variant=="complete"?block.matrix[i][j]-1600.*mass[i][j]:
                (variant=="diffusion"?diffusion[i][j]:(variant=="destruction"?reaction[i][j]:reaction[i][j]+diffusion[i][j]));
            a[i][j]=mass[i][j]/h+l;
        }
        if(variant!="complete")for(int i=0;i<n;++i)
        {
            r[i]=0;const int id=block.rows[i];
            if(variant!="diffusion")r[i]+=resolution.terms.at("inner_k_destruction")[id]+resolution.terms.at("outer_k_destruction")[id];
            if(variant!="destruction")r[i]+=resolution.terms.at("k_diffusion")[id];
        }
        Vector rhs=r;for(auto &v:rhs)v=-v;++solves;const auto delta=Factor(a).Apply(rhs);
        auto linear=Multiply(a,delta);for(int i=0;i<n;++i)linear[i]+=r[i];
        auto raw=s;for(int i=0;i<n;++i)raw[block.rows[i]]+=delta[i];
        const auto bounded=ApplyPimpleKBound(raw,s,time,data);Vector correction(n),rawRate(n),netRate(n),error(n),fullCorrection(size[K()]);
        for(int i=0;i<n;++i){correction[i]=bounded.state[block.rows[i]]-raw[block.rows[i]];rawRate[i]=delta[i]/h;netRate[i]=rawRate[i]+correction[i]/h;error[i]=rawRate[i]-nativeRate[i];fullCorrection[block.rows[i]-offset[K()]]=correction[i];}
        auto storage=Multiply(mass,correction),added=Multiply(a,correction);for(auto &v:storage)v/=h;
        Vector corrected=linear;for(int i=0;i<n;++i)corrected[i]+=added[i];
        std::vector<std::string> samples;
        Vector fullRaw(size[K()]),fullNet(size[K()]);for(int i=0;i<n;++i){fullRaw[block.rows[i]-offset[K()]]=rawRate[i];fullNet[block.rows[i]-offset[K()]]=netRate[i];}
        double positive=0,negative=0;for(const auto &q:Quadrature(K())){const double v=q.weight*q.Value(fullCorrection)/h;(v>=0?positive:negative)+=v;}
        for(const auto &p:points){const auto q=Point(K(),p);samples.push_back(J::Object({{"coordinate",J::Numbers(p)},
            {"raw_rate",J::Number(q.Value(fullRaw))},{"bound_rate",J::Number(q.Value(fullCorrection)/h)},
            {"net_rate",J::Number(q.Value(fullNet))},{"independent_full_spatial_rate",J::Number(q.Value(rate.rate))}}));}
        scans.push_back(J::Object({{"variant",J::Quote(variant)},{"h",J::Number(h)},{"h_sampled_max_destruction_rate",J::Number(h*maximumRate)},
            {"raw_rate_mass_norm",J::Number(massNorm(rawRate))},{"bound_rate_mass_norm",J::Number(allNorm(fullCorrection)/h)},
            {"net_rate_mass_norm",J::Number(massNorm(netRate))},{"raw_minus_full_spatial_rate_mass_norm",J::Number(massNorm(error))},
            {"raw_linear_residual",J::Number(Maximum(linear))},{"corrected_frozen_residual",J::Number(Maximum(corrected))},
            {"storage_reaction",J::Numbers(storage)},{"added_frozen_residual",J::Numbers(added)},
            {"signed_integral_per_time",J::Number(integral(fullCorrection,integralLoad)/h)},
            {"near_integral_per_time",J::Number(integral(fullCorrection,nearLoad)/h)},
            {"positive_negative_integral_per_time",J::Numbers(Vector{positive,negative})},
            {"points",J::List(samples)},{"bound",bounded.report}}));
    }
    const auto sample=PhysicalSamples(s,time,data);double volume=0;for(const auto &q:Quadrature(K()))volume+=q.weight;
    Require(hash==RANSWallR3::Hash(s),"contact diagnostic modified retained state");
    std::vector<std::pair<std::string,std::string>> result{{"schema","1"},{"operation",J::Quote("contact")},{"time",J::Number(time)},
        {"state_sha256",J::Quote(hash)},{"retained_state_unchanged","true"},{"alpha","0"},
        {"spatial_block_norms",J::Numbers(Norms(rawSpatial,constraints,m_robustUstar))},
        {"epsilon_trace_defect",J::Number(trace)},{"contacts",J::List(contacts)},
        {"component_integral_rates",J::List(componentIntegrals)},{"physical_production_destruction_global_near",J::Numbers(productionDestruction)},
        {"physical_outward_advection_diffusion_global",J::Numbers(flux(upper[2]))},{"physical_outward_advection_diffusion_near",J::Numbers(flux(nearHeight))},
        {"near_height",J::Number(nearHeight)},{"global_k_integral",J::Number(integral(kval,integralLoad))},{"near_k_integral",J::Number(integral(kval,nearLoad))},
        {"field_L2_norms_u_v_k_epsilon_pinner_pouter",J::Numbers(fieldNorms)},{"field_integrals",J::Numbers(fieldIntegrals)},
        {"bulk_u",J::Number(fieldIntegrals[0]/volume)},{"flow_through_time",J::Number((upper[0]-lower[0])*volume/fieldIntegrals[0])},
        {"H_over_ustar",J::Number(m_parameters.geometry.height/m_robustUstar)},{"profiles",J::List(profiles)},
        {"wall_stress",J::Numbers(sample.at("wall_stress"))},{"wall_B",J::Numbers(sample.at("wall_B"))},{"wall_lambda",J::Numbers(sample.at("wall_lambda"))},
        {"maximum_sampled_destruction_rate",J::Number(maximumRate)},{"spatial_rate_mass_norm",J::Number(allNorm(rate.rate))},
        {"mass_solve_residual",J::Number(Maximum(rate.freeMassResidual))},{"native_frozen_load_agreement",J::Number(agreement)},
        {"scan",J::List(scans)},{"diagnostic_k_solves",J::Number(solves)},{"production_governing_solves","0"},
        {"scope",J::Quote("fixed state and evaluation time; k history equals this k; no flow/epsilon update, no Lift of selected state; flux estimates need not close a local CG control-volume budget")}};
    if(exportMass){result.push_back({"free_k_mass",J::Rows(mass)});result.push_back({"free_k_rows",J::Numbers(block.rows)});result.push_back({"offsets",J::Numbers(offset)});result.push_back({"k_integral_load",J::Numbers(integralLoad)});result.push_back({"near_k_integral_load",J::Numbers(nearLoad)});}
    return J::Object(result);
}
}
