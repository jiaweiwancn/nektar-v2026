// SPDX-License-Identifier: MIT
// Unforced flat equilibrium qualification. Historical prescribed-top diagnostics
// retain their guards and names; no projected load is used as physical forcing.
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionResolution.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <cmath>
#include <numeric>

namespace Nektar
{
using namespace SolverUtils;
using namespace RANSWallAlgebra;
namespace J=RANSWallSolveJson;

std::string RANSWallEvolutionPatch::EquilibriumSpatialChecks(const Vector &s,
    const RANSWallEvolutionData &data,int order) const
{
    Require(EquilibriumOpen() && !Mapped() && order>=16 && order<=48,"WP2 spatial scope/order");
    const double H=m_parameters.geometry.height,us=m_robustUstar,K0=InflowAtHeight(H).k,tau=us*us;
    const auto wall=Profile(K0,{lower[0],0,0},1.);const auto pi=EquilibriumReferencePressures();
    const Vector bandEdges{0,H,.012,.025,.05,.1,.2,.4,upper[2]};
    const auto parts=ResidualResolutionComponents(s,s,0,0,data,false,false,false,true);
    const auto exact=ResidualResolutionComponents(s,s,0,0,data,false,false,true,true);
    auto r=parts.terms.at("spatial_sum"),ar=exact.terms.at("spatial_sum");
    const auto native=Residual(s,s,0,0,data,false,true);double nativeError=0;
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
    // Columns below carry physical units; full-gradient norms retain all x/z
    // entries, including gradients whose analytical reference is zero.
    const std::vector<std::string> names{"velocity","W","k","epsilon","nut","grad_velocity","grad_k","grad_epsilon","pi","mechanical_pressure","shear","normal_stress","production","divergence","U_x","W_x","W_z"};
    std::array<Vector,2> error{Vector(names.size()),Vector(names.size())},reference=error,maximum=error;
    Vector volume(2),divMaximum(2),extraCovariance(2,1e300),extraScalarMinimum(2,1e300);
    std::vector<std::string> covarianceWitness(2,"null");
    std::map<std::pair<double,double>,Vector> radial,bands;
    for(const auto &q:Quadrature(0,order))
    {
        const double z=q.coordinate[2];const bool inner=z<H;const int region=inner?0:1;
        const auto fixed=InflowAtHeight(std::max(H,z));const double U=EquilibriumReferenceValue(0,q.coordinate);
        const double kr=EquilibriumReferenceValue(K(),q.coordinate),er=inner?wall->Dissipation(z):fixed.epsilon;
        const double du=inner?wall->Original().VelocityGradient(z,{InflowAtHeight(H).U,0},{0,0})[0]:fixed.dU;
        const double dkr=inner?wall->ReferenceKGradient(z):0.,der=inner?0.:-er/z;
        RANSTensor g{},gr{};g[0]=q.Gradient(s,offset[0]);g[2]=q.Gradient(s,offset[1]);gr[0][2]=du;
        const auto gk=q.Gradient(s,offset[K()]);const double k=q.Value(s,offset[K()]);
        const auto eRow=inner?RANSWallGlobalPoint{}:GeometryPoint(Epsilon(),q.geometry,q.coordinate);
        const double epsilon=inner?StateProfile(s,q.coordinate,0,data)->Dissipation(z):eRow.Value(s,offset[Epsilon()]);
        const auto ge=inner?RANSVector{}:eRow.Gradient(s,offset[Epsilon()]);
        RANSTensor T{},Tr{},R{};double nut=0,nutr=0,P=0,Pr=0;
        const auto pRow=GeometryPoint(inner?InnerPressure():OuterPressure(),q.geometry,q.coordinate);
        const double pressure=pRow.Value(s,offset[inner?InnerPressure():OuterPressure()]);
        if(inner)
        {
            nut=StateProfile(s,q.coordinate,0,data)->Original().EddyViscosity(z);nutr=wall->Original().EddyViscosity(z);
            const auto c=InnerStress(k,nut,g),cr=InnerStress(kr,nutr,gr);P=c.production;Pr=cr.production;
            R=c.reynoldsStress;T=RANSWallInnerMeanStress(c,m_parameters.geometry.molecularNu,pressure);Tr=RANSWallInnerMeanStress(cr,m_parameters.geometry.molecularNu,pi[0]);
        }
        else
        {
            const auto c=OuterClosure(k,epsilon,g,q.coordinate),cr=OuterClosure(kr,er,gr,q.coordinate);
            R=c.reynoldsStress;nut=c.eddyViscosity;nutr=cr.eddyViscosity;P=c.production;Pr=cr.production;
            for(int a=0;a<3;++a)for(int b=0;b<3;++b)
            {T[a][b]=-c.reynoldsStress[a][b]+(a==b?2*k/3-pressure:0.);Tr[a][b]=-cr.reynoldsStress[a][b]+(a==b?2*kr/3-pi[1]:0.);}
        }
        const auto eigen=J::Eigenvalues(R);
        if(eigen[0]<extraCovariance[region])
        {extraCovariance[region]=eigen[0];covarianceWitness[region]=J::Object({{"coordinate",J::Numbers(q.coordinate)},{"geometry",J::Number(q.geometry)},
            {"element",J::Number(q.element)},{"region",J::Number(region)},{"radial_side",J::Number(inner?-1:1)},
            {"gradient",J::Rows(g)},{"k",J::Number(k)},{"epsilon",J::Number(epsilon)},{"tensor",J::Rows(R)},{"eigenvalues",J::Numbers(eigen)}});}
        extraScalarMinimum[region]=std::min(extraScalarMinimum[region],k);Require(k>0 && epsilon>0 && std::isfinite(k) && std::isfinite(epsilon),"independent WP2 scalar failure");
        const double u=q.Value(s,offset[0]),w=q.Value(s,offset[1]);double dg=0,gr2=0,kg=0,eg=0;
        for(int a=0;a<3;++a)for(int b=0;b<3;++b){dg+=std::pow(g[a][b]-gr[a][b],2);gr2+=gr[a][b]*gr[a][b];}
        for(int a=0;a<3;++a){kg+=std::pow(gk[a]-(a==2?dkr:0.),2);eg+=std::pow(ge[a]-(a==2?der:0.),2);}
        const double div=g[0][0]+g[2][2];divMaximum[region]=std::max(divMaximum[region],std::abs(div));
        const Vector d{std::hypot(u-U,w),w,k-kr,inner?0:epsilon-er,nut-nutr,std::sqrt(dg),std::sqrt(kg),std::sqrt(eg),pressure-pi[region],pressure-2*k/3-(pi[region]-2*kr/3),T[0][2]-Tr[0][2],T[2][2]-Tr[2][2],P-Pr,div,g[0][0],g[2][0],g[2][2]};
        const Vector v{U,us,kr,inner?0:er,nutr,std::sqrt(gr2),std::abs(dkr),std::abs(der),tau,tau,tau,tau,std::max(std::abs(Pr),1e-30),us/H,us/H,us/H,us/H};
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
    for(const auto &name:RANSWallResolutionSpatialTerms())
    {auto v=parts.terms.at(name),av=exact.terms.at(name);for(auto [i,value]:parts.constraints){(void)value;v[i]=av[i]=0.;}
     components.push_back({name,J::Numbers(Norms(v,parts.constraints,us))});analyticComponents.push_back({name,J::Numbers(Norms(av,parts.constraints,us))});}
    for(int b=0;b<Gauge();++b)
    {std::vector<std::string> elements;for(const auto &e:fields[b]->Elements())elements.push_back(J::Object({{"geometry",J::Number(e.geometry)},{"degree",J::Number(e.degree)},{"lower",J::Numbers(e.lower)},{"upper",J::Numbers(e.upper)}}));mesh.push_back(J::Object({{"block",J::Number(b)},{"size",J::Number(size[b])},{"elements",J::List(elements)}}));}
    double wallError=0,wallMeasure=0,wallMean=0;
    for(const auto &f:SurfaceView(0,2,0.,-1,order))
    {RANSTensor g{};g[0]=f.test.Gradient(s,offset[0]);g[2]=f.test.Gradient(s,offset[1]);
     const double nt=StateProfile(s,f.test.coordinate,0,data)->Original().EddyViscosity(0.);
     const double traction=(m_parameters.geometry.molecularNu+nt)*(g[0][2]+g[2][0]);
     wallError+=f.weight*std::pow(traction-tau,2);wallMean+=f.weight*traction;wallMeasure+=f.weight;}
    const auto audit=CompleteRestrictedAudit(s,s,0,0,data);
    const auto domain=ScalarDomain(s,0,data,true);
    const double divRMS=std::sqrt((error[0][13]+error[1][13])/(volume[0]+volume[1]));
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
        {"wall_shear_RMS_relative_error",J::Number(std::sqrt(wallError/wallMeasure)/tau)},
        {"wall_mean_shear",J::Number(wallMean/wallMeasure)},{"strong_divergence_RMS_over_ustar_H",J::Number(divRMS/(us/H))},
        {"scalar_domain",ScalarDomainJSON(domain,true)},
        {"covariance_inner_outer",J::Numbers(audit.at("accepted_awf_covariance"))},
        {"interface_raw",J::Numbers(audit.at("interface"))},{"pressure_reference_inner_outer",J::Numbers(pi)},
        {"boundary_status",DevelopingChecks(s,0,data)},{"governing_solves","0"},{"physical_steps","0"},
        {"scope",J::Quote("full deployable free-row physical mass-dual residual; analytical loads use analytical boundary values; complete gradients; no stationary correction, anchoring or adoption")}});
}
}
