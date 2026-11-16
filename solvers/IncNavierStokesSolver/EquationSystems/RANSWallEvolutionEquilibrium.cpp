// SPDX-License-Identifier: MIT
// Unforced flat equilibrium qualification. Historical prescribed-top diagnostics
// retain their guards and names; no projected load is used as physical forcing.
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionResolution.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallTiming.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <cmath>
#include <numeric>

namespace Nektar
{
using namespace SolverUtils;
using namespace RANSWallAlgebra;
namespace J=RANSWallSolveJson;

std::string RANSWallEvolutionPatch::EquilibriumSpatialChecks(const Vector &s,
    const RANSWallEvolutionData &data,int order,bool monitoring) const
{
    RANSWallTiming::Scope timing(monitoring?"output.wp2_per_state_monitoring":"diagnostic.wp2_spatial");
    Require(EquilibriumOpen() && !Mapped() && order>=16 && order<=48,"WP2 spatial scope/order");
    const double H=m_parameters.geometry.height,us=m_robustUstar,K0=InflowAtHeight(H).k,tau=us*us;
    const auto wall=Profile(K0,{lower[0],0,0},1.);const auto pi=EquilibriumReferencePressures();
    const Vector bandEdges{0,H,.012,.025,.05,.1,.2,.4,upper[2]};
    ResidualResolutionResult parts,exact;
    if(!monitoring)
    {parts=ResidualResolutionComponents(s,s,0,0,data,false,false,false,true);
     exact=ResidualResolutionComponents(s,s,0,0,data,false,false,true,true);}
    else
    {parts.constraints=Constraints(s,0,data);parts.terms["spatial_sum"]=Residual(s,s,0,0,data,false,true);
     exact.terms["spatial_sum"]=Vector(Size());}
    auto r=parts.terms.at("spatial_sum"),ar=exact.terms.at("spatial_sum");
    const auto native=monitoring?r:Residual(s,s,0,0,data,false,true);double nativeError=0;
    std::vector<std::string> freeRows;Vector maxRow(4);
    for(int b=0;b<InnerPressure();++b)
    {
        Vector rows;const double scale=b<2?us:b==K()?K0:std::pow(us,3)/(m_parameters.model.Kappa()*H);
        const double rate=scale*us/H;
        for(int i=0;i<size[b];++i)if(!parts.constraints.count(offset[b]+i))
        {rows.push_back(offset[b]+i);nativeError=std::max(nativeError,std::abs(r[offset[b]+i]-native[offset[b]+i]));
         maxRow[b]=std::max(maxRow[b],std::abs(r[offset[b]+i])/(rate*std::sqrt(m_volumes[b]*m_mass[b][i][i])));}
        freeRows.push_back(J::Numbers(rows));
    }
    for(auto [i,value]:parts.constraints){(void)value;r[i]=ar[i]=0.;}
    const auto norms=Norms(r,parts.constraints,us),anorms=Norms(ar,parts.constraints,us);
    // Local loads still share global CG rows. Contract each panel against the
    // inverse of the assembled free-row mass and total load, retaining every
    // cross term. These signed projections sum to the total norm; they are not
    // positive error allocations or percentages of a traction discrepancy.
    std::vector<std::pair<std::string,std::string>> signedPanels;
    Vector projectionSum(4);
    if(!monitoring)
    {
    std::map<std::string,Vector> panelLoads;
    for(const auto &[name,load]:parts.terms)if(name.rfind("radial:",0)==0)
    {
        const auto key=name.substr(0,name.find(':',name.find(':',7)+1));
        auto &v=panelLoads[key];if(v.empty())v=Vector(Size());
        for(int i=0;i<Size();++i)v[i]+=load[i];
    }
    Vector remaining=r;
    for(const auto &[name,load]:panelLoads)for(int i=0;i<Size();++i)remaining[i]-=load[i];
    panelLoads["boundary_interface_and_other_loads"]=remaining;
    std::map<std::string,Vector> panelProjections;
    for(int b=0;b<InnerPressure();++b)
    {
        std::vector<int> free;Vector total;
        for(int i=0;i<size[b];++i)if(!parts.constraints.count(offset[b]+i))
        {free.push_back(i);total.push_back(r[offset[b]+i]);}
        WallMatrix mass(free.size(),Vector(free.size()));
        for(int i=0;i<free.size();++i)for(int j=0;j<free.size();++j)mass[i][j]=m_mass[b][free[i]][free[j]];
        const Factor factor(mass);const auto dual=factor.Apply(total);
        const double dimNorm=std::sqrt(std::max(0.,std::inner_product(total.begin(),total.end(),dual.begin(),0.)));
        const double scale=b<2?us:b==K()?K0:std::pow(us,3)/(m_parameters.model.Kappa()*H);
        const double normalization=(scale*us/H)*std::sqrt(m_volumes[b]);
        for(const auto &[name,load]:panelLoads)
        {
            long double dot=0;for(int j=0;j<free.size();++j)dot+=(long double)load[offset[b]+free[j]]*dual[j];
            auto &v=panelProjections[name];if(v.empty())v=Vector(4);
            v[b]=dimNorm>1e-30?double(dot)/(dimNorm*normalization):0.;projectionSum[b]+=v[b];
        }
        Require(std::abs(projectionSum[b]-norms[b])<1e-10*std::max(1.,norms[b]),"signed radial projections do not close assembled residual");
    }
    for(const auto &[name,value]:panelProjections)signedPanels.push_back({name,J::Numbers(value)});
    }
    // Columns below carry physical units; full-gradient norms retain all x/z
    // entries, including gradients whose analytical reference is zero.
    const std::vector<std::string> names{"velocity","W","k","epsilon","nut","grad_velocity","grad_k","grad_epsilon","pi","mechanical_pressure","shear","normal_stress","production","divergence","U_x","W_x","W_z"};
    std::array<Vector,2> error{Vector(names.size()),Vector(names.size())},reference=error,maximum=error;
    Vector volume(2),divMaximum(2),extraCovariance(2,1e300),extraScalarMinimum(2,1e300);
    std::vector<std::string> covarianceWitness(2,"null");
    std::map<std::pair<double,double>,Vector> radial,bands;
    auto observations=Quadrature(0,order);
    // Zero integration weights retain each physical station and both native
    // radial/tangential limits without changing any physical norm or load.
    const auto rule=LibUtilities::PointsManager()[LibUtilities::PointsKey(order,LibUtilities::eGaussGaussLegendre)];
    for(const auto &el:fields[0]->Elements())for(double x:Vector{-.2,-.15,-.1,-.05,0.,.05,.1,.3,.6})
        if(x>=el.lower[0] && x<=el.upper[0])for(int j=-1;j<=order;++j)
        {
            const double z=j<0?el.lower[2]:j==order?el.upper[2]:el.lower[2]+(rule->GetZ()[j]+1)*(el.upper[2]-el.lower[2])/2;
            auto point=GeometryPoint(0,el.geometry,{x,0,z});point.weight=0;observations.push_back(std::move(point));
        }
    for(const auto &q:observations)
    {
        const double z=q.coordinate[2];const bool inner=fields[0]->Elements()[q.element].upper[2]<=H;const int region=inner?0:1;
        const auto &ref=EquilibriumReferenceSample(inner,z);
        const double U=ref[0],kr=ref[1],er=ref[2],du=ref[3],dkr=ref[4],der=ref[5];
        RANSTensor g{},gr{};g[0]=q.Gradient(s,offset[0]);g[2]=q.Gradient(s,offset[1]);gr[0][2]=du;
        const auto gk=q.Gradient(s,offset[K()]);const double k=q.Value(s,offset[K()]);
        const auto eRow=inner?RANSWallGlobalPoint{}:GeometryPoint(Epsilon(),q.geometry,q.coordinate);
        const double epsilon=inner?StateProfile(s,q.coordinate,0,data)->Dissipation(z):eRow.Value(s,offset[Epsilon()]);
        const auto ge=inner?RANSVector{}:eRow.Gradient(s,offset[Epsilon()]);
        RANSTensor T{},Tr{},R{};Tr[0][0]=ref[8];Tr[2][2]=ref[9];Tr[0][2]=ref[10];
        double nut=0,P=0;const double nutr=ref[6],Pr=ref[7];
        const auto pRow=GeometryPoint(inner?InnerPressure():OuterPressure(),q.geometry,q.coordinate);
        const double pressure=pRow.Value(s,offset[inner?InnerPressure():OuterPressure()]);
        if(inner)
        {
            nut=StateProfile(s,q.coordinate,0,data)->Original().EddyViscosity(z);
            const auto c=InnerStress(k,nut,g);P=c.production;
            R=c.reynoldsStress;T=RANSWallInnerMeanStress(c,m_parameters.geometry.molecularNu,pressure);
        }
        else
        {
            const auto c=OuterClosure(k,epsilon,g,q.coordinate);
            R=c.reynoldsStress;nut=c.eddyViscosity;P=c.production;
            for(int a=0;a<3;++a)for(int b=0;b<3;++b)
            T[a][b]=-c.reynoldsStress[a][b]+(a==b?2*k/3-pressure:0.);
        }
        const auto eigen=J::Eigenvalues(R);
        if(eigen[0]<extraCovariance[region])
        {extraCovariance[region]=eigen[0];covarianceWitness[region]=J::Object({{"coordinate",J::Numbers(q.coordinate)},{"geometry",J::Number(q.geometry)},
            {"element",J::Number(q.element)},{"region",J::Number(region)},{"radial_side",J::Number(z==H?(inner?-1:1):0)},
            {"epsilon_context",J::Quote(inner?"dependent inner E_C at transported endpoint K":"actual transported outer epsilon")},
            {"gradient",J::Rows(g)},{"k",J::Number(k)},{"epsilon",J::Number(epsilon)},{"tensor",J::Rows(R)},{"eigenvalues",J::Numbers(eigen)}});}
        extraScalarMinimum[region]=std::min(extraScalarMinimum[region],k);Require(k>0 && epsilon>0 && std::isfinite(k) && std::isfinite(epsilon),"independent WP2 scalar failure");
        const double u=q.Value(s,offset[0]),w=q.Value(s,offset[1]);double dg=0,gr2=0,kg=0,eg=0;
        for(int a=0;a<3;++a)for(int b=0;b<3;++b){dg+=std::pow(g[a][b]-gr[a][b],2);gr2+=gr[a][b]*gr[a][b];}
        for(int a=0;a<3;++a){kg+=std::pow(gk[a]-(a==2?dkr:0.),2);eg+=std::pow(ge[a]-(a==2?der:0.),2);}
        const double div=g[0][0]+g[2][2];divMaximum[region]=std::max(divMaximum[region],std::abs(div));
        const Vector d{std::hypot(u-U,w),w,k-kr,inner?0:epsilon-er,nut-nutr,std::sqrt(dg),std::sqrt(kg),std::sqrt(eg),pressure-pi[region],pressure-2*k/3-(pi[region]-2*kr/3),T[0][2]-Tr[0][2],T[2][2]-Tr[2][2],P-Pr,div,g[0][0],g[2][0],g[2][2]};
        const Vector v{U,us,kr,inner?0:er,nutr,std::sqrt(gr2),inner?std::abs(dkr):K0/H,inner?0:std::abs(der),tau,tau,tau,tau,std::max(std::abs(Pr),1e-30),us/H,us/H,us/H,us/H};
        const auto &el=fields[0]->Elements()[q.element];auto &panel=radial[{el.lower[2],el.upper[2]}];if(panel.empty())panel=Vector(names.size());
        volume[region]+=q.weight;
        for(int i=0;i<names.size();++i)
        {error[region][i]+=q.weight*d[i]*d[i];reference[region][i]+=q.weight*v[i]*v[i];panel[i]+=q.weight*d[i]*d[i];
         if(i!=0 && (i!=3 || !inner))maximum[region][i]=std::max(maximum[region][i],std::abs(d[i])/std::max(std::abs(v[i]),1e-30));}
        maximum[region][0]=std::max(maximum[region][0],std::abs(u-U)/std::max(std::abs(U),us));
        for(int j=1;j<bandEdges.size();++j)if(z>=bandEdges[j-1] && z<=bandEdges[j])
        {auto &v=bands[{bandEdges[j-1],bandEdges[j]}];if(v.empty())v=Vector(6);v[0]+=q.weight;
         v[1]+=q.weight*std::pow((u-U)/std::max(std::abs(U),us),2);v[2]+=q.weight*std::pow((k-kr)/kr,2);
         v[3]+=q.weight*(inner?0.:std::pow((epsilon-er)/er,2));v[4]+=q.weight*std::pow((nut-nutr)/std::max(nutr,1e-30),2);v[5]+=q.weight*std::pow(w/us,2);break;}
    }
    std::vector<std::string> regions,panels,mesh,bandErrors;
    std::vector<std::pair<std::string,std::string>> components,analyticComponents;
    for(int region=0;region<2;++region)
    {Vector raw=error[region],relative=raw;for(int i=0;i<raw.size();++i){raw[i]=std::sqrt(raw[i]);relative[i]=reference[region][i]>1e-60?raw[i]/std::sqrt(reference[region][i]):0.;}
     regions.push_back(J::Object({{"volume",J::Number(volume[region])},{"dimensional_L2",J::Numbers(raw)},{"relative_L2",J::Numbers(relative)},
        {"pointwise_normalized_maxima",J::Numbers(maximum[region])},{"strong_divergence_max",J::Number(divMaximum[region])},
        {"zero_reference_gradient_scale",J::Quote("report dimensional grad-k and streamwise gradient norms; u_tau/H for zero velocity-gradient references; inner epsilon is dependent, not transported")}}));}
    for(auto [bounds,v]:radial){for(auto &x:v)x=std::sqrt(x);panels.push_back(J::Object({{"bounds",J::Numbers(Vector{bounds.first,bounds.second})},{"dimensional_L2",J::Numbers(v)}}));}
    for(auto [bounds,v]:bands){const double measure=v[0];for(int j=1;j<v.size();++j)v[j]=std::sqrt(v[j]/measure);
        bandErrors.push_back(J::Object({{"bounds",J::Numbers(Vector{bounds.first,bounds.second})},{"volume_then_local_relative_RMS_U_k_epsilon_nut_W",J::Numbers(v)}}));}
    if(!monitoring)for(const auto &name:RANSWallResolutionSpatialTerms())
    {auto v=parts.terms.at(name),av=exact.terms.at(name);for(auto [i,value]:parts.constraints){(void)value;v[i]=av[i]=0.;}
     components.push_back({name,J::Numbers(Norms(v,parts.constraints,us))});analyticComponents.push_back({name,J::Numbers(Norms(av,parts.constraints,us))});}
    if(!monitoring)for(int b=0;b<Gauge();++b)
    {std::vector<std::string> elements;for(const auto &e:fields[b]->Elements())elements.push_back(J::Object({{"geometry",J::Number(e.geometry)},{"degree",J::Number(e.degree)},{"lower",J::Numbers(e.lower)},{"upper",J::Numbers(e.upper)}}));mesh.push_back(J::Object({{"block",J::Number(b)},{"size",J::Number(size[b])},{"elements",J::List(elements)}}));}
    double wallError=0,wallMeasure=0,wallMean=0;
    for(const auto &f:SurfaceView(0,2,0.,-1,order))
    {RANSTensor g{};g[0]=f.test.Gradient(s,offset[0]);g[2]=f.test.Gradient(s,offset[1]);
     const double nt=StateProfile(s,f.test.coordinate,0,data)->Original().EddyViscosity(0.);
     const double traction=(m_parameters.geometry.molecularNu+nt)*(g[0][2]+g[2][0]);
     wallError+=f.weight*std::pow(traction-tau,2);wallMean+=f.weight*traction;wallMeasure+=f.weight;}
    const double divRMS=std::sqrt((error[0][13]+error[1][13])/(volume[0]+volume[1]));
    if(monitoring)return J::Object({{"version",J::Quote("WP2PerAcceptedStateV1")},
        {"state_sha256",J::Quote(RANSWallR3::Hash(s))},{"unanchored_free_block_norms",J::Numbers(norms)},
        {"maximum_scaled_individual_row",J::Numbers(maxRow)},{"regional_errors_inner_outer",J::List(regions)},
        {"log_height_bands",J::List(bandErrors)},{"volume_covariance_inner_outer",J::Numbers(extraCovariance)},
        {"volume_minimum_k_inner_outer",J::Numbers(extraScalarMinimum)},
        {"wall_shear_RMS_relative_error",J::Number(std::sqrt(wallError/wallMeasure)/tau)},
        {"wall_mean_shear",J::Number(wallMean/wallMeasure)},{"strong_divergence_RMS_over_ustar_H",J::Number(divRMS/(us/H))},
        {"pressure_reference_inner_outer",J::Numbers(pi)},{"boundary_status",DevelopingChecks(s,0,data)},
        {"reference_cache_height_entries",J::Number(m_equilibriumReferences.size())},
        {"scope",J::Quote("every retained accepted state; physical GL16 on all native panels; complete production sided audits remain in state_audits; fixed analytical reference cache only")}});
    const auto audit=CompleteRestrictedAudit(s,s,0,0,data);
    const auto domain=ScalarDomain(s,0,data,true);
    return J::Object({{"version",J::Quote("WP2UnforcedSpatialV1")},{"state_sha256",J::Quote(RANSWallR3::Hash(s))},{"state",J::Numbers(s)},
        {"size",J::Number(Size())},{"sizes",J::Numbers(size)},{"mesh",J::List(mesh)},{"quadrature_order",J::Number(order)},
        {"free_coefficient_rows",J::List(freeRows)},{"unanchored_free_block_norms",J::Numbers(norms)},
        {"analytical_free_block_norms",J::Numbers(anorms)},{"maximum_scaled_individual_row",J::Numbers(maxRow)},
        {"native_decomposition_free_max_error",J::Number(nativeError)},{"metric_columns",J::Quote("velocity,W,k,epsilon,nut,grad_velocity,grad_k,grad_epsilon,pi,mechanical_pressure,shear,normal_stress,production,divergence,U_x,W_x,W_z")},
        {"regional_errors_inner_outer",J::List(regions)},{"radial_errors",J::List(panels)},
        {"log_height_bands",J::List(bandErrors)},
        {"independent_covariance_inner_outer",J::Numbers(extraCovariance)},{"independent_covariance_witnesses",J::List(covarianceWitness)},
        {"independent_minimum_k_inner_outer",J::Numbers(extraScalarMinimum)},
        {"global_components",J::Object(components)},{"analytical_components",J::Object(analyticComponents)},
        {"signed_radial_mass_dual_projections",J::Object(signedPanels)},{"signed_projection_sum",J::Numbers(projectionSum)},
        {"signed_projection_definition",J::Quote("panel_load^T M_free^-1 total_load / norm(total_load) / native rate-volume scale; shared rows and cross terms retained")},
        {"wall_shear_RMS_relative_error",J::Number(std::sqrt(wallError/wallMeasure)/tau)},
        {"wall_mean_shear",J::Number(wallMean/wallMeasure)},{"strong_divergence_RMS_over_ustar_H",J::Number(divRMS/(us/H))},
        {"scalar_domain",ScalarDomainJSON(domain,true)},
        {"covariance_inner_outer",J::Numbers(audit.at("accepted_awf_covariance"))},
        {"interface_raw",J::Numbers(audit.at("interface"))},{"pressure_reference_inner_outer",J::Numbers(pi)},
        {"boundary_status",DevelopingChecks(s,0,data)},{"governing_solves","0"},{"physical_steps","0"},
        {"scope",J::Quote("full deployable free-row physical mass-dual residual; analytical loads use analytical boundary values; complete gradients; no stationary correction, anchoring or adoption")}});
}

std::string RANSWallEvolutionPatch::ResourcePolicyJSON() const
{
    return J::Object({{"version",J::Quote(m_resourcePolicy)},{"dense_slot_limit",J::Number(m_denseSizeLimit)},
        {"mass",J::Quote("unchanged physical consistent mass; B=(alpha*M_free)^-1")},
        {"scope",J::Quote(m_resourcePolicy=="Dense3000V1"?"historical production default":"explicit measured WP2 affine serial one-sweep case only; 24GiB RSS / aggregate wall envelope")}});
}

const RANSWallEvolutionPatch::Vector &RANSWallEvolutionPatch::EquilibriumReferenceSample(bool inner,double z) const
{
    Require(EquilibriumOpen() && !Mapped(),"immutable analytical sampling scope");
    const int region=inner?0:1;const auto key=std::make_pair(region,z);const auto old=m_equilibriumReferences.find(key);
    if(old!=m_equilibriumReferences.end())return old->second;
    const double H=m_parameters.geometry.height;const auto pi=EquilibriumReferencePressures();const auto fixed=InflowAtHeight(std::max(H,z));
    const auto wall=Profile(InflowAtHeight(H).k,{lower[0],0,0},1.);
    const double U=inner?wall->Original().Velocity(z,{InflowAtHeight(H).U,0},{0,0})[0]:fixed.U;
    const double k=inner?wall->ReferenceK(z):fixed.k,e=inner?wall->Dissipation(z):fixed.epsilon;
    const double du=inner?wall->Original().VelocityGradient(z,{InflowAtHeight(H).U,0},{0,0})[0]:fixed.dU;
    const double dk=inner?wall->ReferenceKGradient(z):0.,de=inner?0.:-e/z;
    RANSTensor g{},T{};g[0][2]=du;double nt,P;
    if(inner)
    {nt=wall->Original().EddyViscosity(z);const auto c=InnerStress(k,nt,g);P=c.production;T=RANSWallInnerMeanStress(c,m_parameters.geometry.molecularNu,pi[0]);}
    else
    {const auto c=CalibratedStandard()?EvaluateRANSClosure(m_parameters.model.Model(),k,e,g,z,m_parameters.model.Parameters(),0):
        EvaluateRANSRobustClosure(m_parameters.model.Model(),k,e,g,z,m_parameters.model.Parameters(),0,RegularizedReferenceAt({lower[0],0,z}));
     nt=c.eddyViscosity;P=c.production;for(int a=0;a<3;++a)for(int b=0;b<3;++b)T[a][b]=-c.reynoldsStress[a][b]+(a==b?2*k/3-pi[1]:0.);}
    return m_equilibriumReferences.emplace(key,Vector{U,k,e,du,dk,de,nt,P,T[0][0],T[2][2],T[0][2]}).first->second;
}

std::string RANSWallEvolutionPatch::EquilibriumCommonSamples(const Vector &s,const RANSWallEvolutionData &data) const
{
    Require(EquilibriumOpen() && !Mapped(),"common physical samples require the WP2 affine case");
    RANSWallTiming::Scope timing("output.wp2_common_samples");
    Vector knots=m_equilibriumReportingZ;
    if(knots.empty())for(const auto &el:fields[0]->Elements()){knots.push_back(el.lower[2]);knots.push_back(el.upper[2]);}
    std::sort(knots.begin(),knots.end());knots.erase(std::unique(knots.begin(),knots.end()),knots.end());
    Require(knots.front()==0 && knots.back()==upper[2],"common reporting partition does not cover the entire physical domain");
    for(const auto &el:fields[0]->Elements())for(double z:{el.lower[2],el.upper[2]})
        Require(std::find(knots.begin(),knots.end(),z)!=knots.end(),"common reporting partition omits a native interface");
    const auto gl=LibUtilities::PointsManager()[LibUtilities::PointsKey(8,LibUtilities::eGaussGaussLegendre)];
    std::vector<std::string> rows,references;Vector weights,regions;std::set<std::pair<double,double>> horizontal;
    for(const auto &el:fields[0]->Elements())horizontal.insert({el.lower[0],el.upper[0]});
    const double H=m_parameters.geometry.height;
    for(const auto &[lo,hi]:horizontal)for(int ix=0;ix<8;++ix)for(int panel=1;panel<knots.size();++panel)for(int iz=0;iz<8;++iz)
    {
        const double x=lo+(gl->GetZ()[ix]+1)*(hi-lo)/2,z=knots[panel-1]+(gl->GetZ()[iz]+1)*(knots[panel]-knots[panel-1])/2;
        const auto q=Point(0,{x,0,z});const bool inner=z<H;const double u=q.Value(s,offset[0]),w=q.Value(s,offset[1]),k=q.Value(s,offset[K()]);
        const auto profile=inner?StateProfile(s,q.coordinate,0,data):nullptr;
        const double e=inner?profile->Dissipation(z):GeometryPoint(Epsilon(),q.geometry,q.coordinate).Value(s,offset[Epsilon()]);
        RANSTensor g{},T{};g[0]=q.Gradient(s,offset[0]);g[2]=q.Gradient(s,offset[1]);
        const int b=inner?InnerPressure():OuterPressure();const double pi=GeometryPoint(b,q.geometry,q.coordinate).Value(s,offset[b]);double nt,P;
        if(inner)
        {nt=profile->Original().EddyViscosity(z);const auto c=InnerStress(k,nt,g);P=c.production;T=RANSWallInnerMeanStress(c,m_parameters.geometry.molecularNu,pi);}
        else
        {const auto c=OuterClosure(k,e,g,q.coordinate);nt=c.eddyViscosity;P=c.production;for(int a=0;a<3;++a)for(int b=0;b<3;++b)T[a][b]=-c.reynoldsStress[a][b]+(a==b?2*k/3-pi:0.);}
        rows.push_back(J::Numbers(Vector{u,w,k,e,pi,pi-2*k/3,nt,T[0][2],T[2][2],P,g[0][0],g[0][2],g[2][0],g[2][2]}));
        const auto &ref=EquilibriumReferenceSample(inner,z);const auto pressure=EquilibriumReferencePressures();const double pr=pressure[inner?0:1];
        references.push_back(J::Numbers(Vector{ref[0],0,ref[1],ref[2],pr,pr-2*ref[1]/3,ref[6],ref[10],ref[9],ref[7],0,ref[3],0,0}));
        weights.push_back(gl->GetW()[ix]*(hi-lo)/2*gl->GetW()[iz]*(knots[panel]-knots[panel-1])/2);regions.push_back(inner?0:1);
    }
    return J::Object({{"version",J::Quote("WP2CommonPhysicalSamplesV1")},{"state_sha256",J::Quote(RANSWallR3::Hash(s))},
        {"physical_z_partition",J::Numbers(knots)},{"order","8"},{"weights",J::Numbers(weights)},{"regions",J::Numbers(regions)},
        {"columns",J::Quote("U,W,k,epsilon,pi,mechanical_pressure,nut,total_shear,total_normal_stress,production,U_x,U_z,W_x,W_z")},
        {"values",J::List(rows)},{"reference_values",J::List(references)},{"pressure_gauge",J::Quote("absolute total-pressure outlet convention; no post-hoc gauge removal")},
        {"epsilon_context",J::Quote("inner dependent E_C; outer actual retained transported coefficients; outer norms use region=1 only")},
        {"scope",J::Quote("same union of native physical panels for coarse/fine/confirmation; cubic field squared differences integrate exactly; rational constitutive differences need quadrature confirmation")}});
}

RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::EquilibriumOuterKPerturbation(
    const Vector &base,const RANSWallEvolutionData &data,std::string &report) const
{
    Require(EquilibriumOpen() && !Mapped(),"WP2 recovery setup scope");
    RANSWallTiming::Scope timing("startup.wp2_perturbation_projection");
    const double zlo=.012,zhi=.05,K0=InflowAtHeight(m_parameters.geometry.height).k,amplitude=.04;
    std::set<int> fixed;
    for(const auto &[row,value]:ExternalConstraints(0,data))if(row>=offset[K()] && row<offset[K()]+size[K()])fixed.insert(row-offset[K()]);
    for(int i=0;i<fields[K()]->Elements().size();++i)
    {
        const auto &el=fields[K()]->Elements()[i];
        if(el.lower[2]<zlo || el.upper[2]>zhi)
        {const auto q=fields[K()]->ElementPoint(i,{(el.lower[0]+el.upper[0])/2,0,(el.lower[2]+el.upper[2])/2});for(int id:q.ids)fixed.insert(id);}
    }
    const auto endRows=fields[K()]->Boundary([](const RANSWallCoordinate&){return 0.;},
        [&](int axis,int side,const RANSWallGlobalField::Element &el){return axis==0 && side==1 && el.upper[0]==upper[0];});
    for(const auto &[id,value]:endRows)fixed.insert(id);
    Vector load(size[K()]);
    const auto shape=[&](const RANSWallCoordinate &p)
    {if(p[2]<=zlo || p[2]>=zhi)return 0.;const double a=std::sin(std::acos(-1.)*(p[0]-lower[0])/(upper[0]-lower[0]));
     const double b=std::sin(std::acos(-1.)*(p[2]-zlo)/(zhi-zlo));return amplitude*K0*a*a*b*b;};
    for(const auto &q:Quadrature(K(),16))if(q.coordinate[2]>=zlo && q.coordinate[2]<=zhi)
        for(int j=0;j<q.ids.size();++j)load[q.ids[j]]+=q.weight*q.phi[j]*shape(q.coordinate);
    std::vector<int> free;Vector rhs;
    for(int i=0;i<size[K()];++i)if(!fixed.count(i)){free.push_back(offset[K()]+i);rhs.push_back(load[i]);}
    const auto delta=RANSWallAlgebra::Solve(PicardKMass(free),rhs);Vector state=base;
    for(int j=0;j<free.size();++j)state[free[j]]+=delta[j];
    std::set<int> allowed(free.begin(),free.end());Vector frozenBase,frozenState;
    for(int i=0;i<Size();++i)if(!allowed.count(i)){frozenBase.push_back(base[i]);frozenState.push_back(state[i]);}
    Require(RANSWallR3::Hash(frozenBase)==RANSWallR3::Hash(frozenState),"WP2 recovery changed protected coefficient bits");
    double square=0,integral=0,maximum=0,trace=0;
    for(const auto &q:Quadrature(K(),24))
    {const double d=q.Value(state,offset[K()])-q.Value(base,offset[K()]);square+=q.weight*d*d;integral+=q.weight*d;maximum=std::max(maximum,std::abs(d)/K0);}
    for(const auto &face:SurfaceView(K(),2,m_parameters.geometry.height,1,16))
        trace=std::max(trace,std::abs(face.test.Value(state,offset[K()])-face.test.Value(base,offset[K()])));
    Require(maximum>=.02 && maximum<.05 && trace==0,"actual recovery perturbation is unresolved, excessive or changes endpoint K");
    report=J::Object({{"mode",J::Quote("WP2OuterKBumpV1")},{"immutable_base_sha256",J::Quote(RANSWallR3::Hash(base))},
        {"seed_sha256",J::Quote(RANSWallR3::Hash(state))},{"reference_defect_boundary",J::Quote("retained unchanged; perturbation is applied after freezing data and before fresh histories; no Lift")},
        {"frozen_rows_sha256",J::Quote(RANSWallR3::Hash(frozenBase))},{"amplitude_relative_K0",J::Number(amplitude)},
        {"support_z",J::Numbers(Vector{zlo,zhi})},{"support_x",J::Numbers(Vector{lower[0],upper[0]})},
        {"actual_increment_L2",J::Number(std::sqrt(square))},{"initial_TKE_increment",J::Number(integral)},
        {"maximum_actual_relative_k_departure",J::Number(maximum)},{"endpoint_K_change",J::Number(trace)},
        {"setup_scalar_mass_projections","1"},{"governing_solves","0"},{"velocity_pressure_epsilon_changes","0"}});
    return state;
}

}
