// SPDX-License-Identifier: MIT
#include <SolverUtils/Turbulence/RANSFiniteInflowReference.h>
#include <SolverUtils/Turbulence/RANSTwoRegionHill.h>
#include <SolverUtils/Turbulence/EquilibriumAWFProfiles.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <fstream>
#include <iostream>
#include <functional>
#include <map>

using namespace Nektar::SolverUtils;
namespace J=Nektar::RANSWallSolveJson;
namespace
{
int checks=0;
void Check(bool b,const char *s){++checks;if(!b)throw std::runtime_error(s);}
template<class F>void Reject(F f,const char *s){bool b=false;try{f();}catch(const std::invalid_argument&){b=true;}Check(b,s);}
double Norm(const RANSVector &v){return std::hypot(v[0],v[1],v[2]);}
RANSVector Sub(const RANSVector&a,const RANSVector&b){return {a[0]-b[0],a[1]-b[1],a[2]-b[2]};}
double TensorNorm(const RANSTensor&t){double a=0;for(auto &r:t)for(double x:r)a+=x*x;return std::sqrt(a);}
std::vector<std::pair<double,double>> Gauss(int n)
{
    std::vector<std::pair<double,double>> q;const double pi=std::acos(-1.);
    for(int i=0;i<n;++i)
    {
        double x=std::cos(pi*(i+.75)/(n+.5)),dp=0;
        for(int j=0;j<30;++j)
        {
            double p0=1,p1=x;for(int k=2;k<=n;++k){double p=((2*k-1)*x*p1-(k-1)*p0)/k;p0=p1;p1=p;}
            dp=n*(x*p1-p0)/(x*x-1);double d=p1/dp;x-=d;if(std::abs(d)<2e-16)break;
        }
        q.emplace_back(x,2/((1-x*x)*dp*dp));
    }
    return q;
}
double Integrate(const std::function<double(double)> &f,double a,double b,int n=64)
{long double v=0;for(auto [x,w]:Gauss(n))v+=w*f((a+b)/2+(b-a)*x/2)*(b-a)/2;return v;}

std::string GeometryChecks()
{
    RANSTwoRegionHill g(.04,.1,.015,.8,-.2,.2),flat(0,.1,.015,.8,-.2,.2);
    Check(std::abs(g.CurvatureBound()-19.73920880217871)<1e-12,"hill curvature");
    Check(std::abs(g.Inner(0,0).point[2]-.04)<1e-15,"hill height");
    const auto footIn=g.Inner(-.1,0,+1),footOut=g.Inner(-.1,0,-1);
    Check(Norm(Sub(footIn.point,footOut.point))==0 && footIn.baseArea==1 && footOut.baseArea==1,"C1 foot position/slope");
    Check(std::abs(footIn.kappa1+footIn.kappa2+g.CurvatureBound())<1e-12 && footOut.kappa1==0 && footOut.kappa2==0,"one-sided foot curvature");
    Reject([]{RANSTwoRegionHill b(.04,.1,.06,.8,-.2,.2);},"invalid reach accepted");
    double inverse=0,jac=0,gradient=0,interface=0,area=0,minJ=1e100,wall=0,wallBound=0;
    for(double s:{-.19,-.1,-.073,-.031,0.,.046,.1,.19})for(double eta:{0.,.23,.8,1.})
    {
        const auto q=g.Outer(s,eta);minJ=std::min(minJ,q.determinant);
        const auto qi=g.OuterInverse(q.point);inverse=std::max(inverse,std::abs(qi[0]-s)+std::abs(qi[2]-eta));
        const auto f=flat.Outer(s,eta);Check(std::abs(f.point[0]-s)<1e-15 && std::abs(f.point[2]-(.015+eta*(.8-.015)))<1e-15,"flat recovery");
        const auto h=g.Inner(s,.015);interface=std::max(interface,Norm(Sub(h.point,g.Outer(s,0).point)));
        const auto ai=g.InterfaceArea(s,true),ao=g.InterfaceArea(s,false);
        area=std::max(area,Norm({ai[0]+ao[0],0,ai[2]+ao[2]}));
        Check(std::abs(Norm(ai)-std::hypot(h.jacobian[0][0],h.jacobian[2][0]))<1e-15,"interface physical measure");
        if(std::abs(std::abs(s)-.1)>1e-12 && eta>.001 && eta<.999)
        {
            const double d=1e-7;const auto a=g.Outer(s+d,eta),b=g.Outer(s-d,eta),c=g.Outer(s,eta+d),e=g.Outer(s,eta-d);
            RANSTensor fd{};for(int t:{0,2}){fd[t][0]=(a.point[t]-b.point[t])/(2*d);fd[t][2]=(c.point[t]-e.point[t])/(2*d);}
            jac=std::max(jac,std::abs(fd[0][0]*fd[2][2]-fd[0][2]*fd[2][0]-q.determinant));
            for(int a:{0,2})
            {
                double grad=0;for(int b:{0,2})grad+=q.inverse[b][a]*(fd[0][b]+2*fd[2][b]);
                gradient=std::max(gradient,std::abs(grad-(a==0?1:2)));
            }
        }
    }
    for(double s:{-.18,-.067,0.,.038,.18})for(double r:{.001,.01,.014})
    {
        const auto p=g.Inner(s,r);const auto inv=g.InnerInverse(p.point);
        inverse=std::max(inverse,Norm(Sub(inv,{s,0,r})));
        const auto distance=g.PhysicalWallDistance(p.point,1e-9);
        wall=std::max(wall,std::abs(distance.distance-r));wallBound=std::max(wallBound,distance.errorBound);
    }
    const auto roofDistance=g.PhysicalWallDistance(g.Outer(.053,.5).point,1e-9);
    Check(std::abs(roofDistance.distance-.5)>1e-3,"outer eta used as physical wall distance");
    double volume=0,volume32=0,flux=0;
    const std::vector<double> panels{-.2,-.1,0,.1,.2};
    for(int i=0;i<4;++i)
    {
        auto integrand=[&](double s)
        {
            const auto q=g.Inner(s,0),h=g.Inner(s,.015);
            return q.baseArea*(.015+(q.kappa1+q.kappa2)*.015*.015/2)+g.Outer(s,.5).determinant;
        };
        volume+=Integrate(integrand,panels[i],panels[i+1]);volume32+=Integrate(integrand,panels[i],panels[i+1],32);
    }
    // Independent physical-graph boundary integral of (x,z): top=.8*width,
    // bottom=int(x*f'-f)dx; sides=xmax*.8-xmin*.8.
    const double bottom=Integrate([](double x){const double a=std::acos(-1.)/.1;const double f=.02*(1+std::cos(a*x)),df=-.02*a*std::sin(a*x);return x*df-f;},-.1,.1);
    flux=.8*.4+.8*.4+bottom;
    const double exact=.4*.8-.04*.1;
    Check(inverse<1e-10 && jac<1e-8 && gradient<1e-8,"map inverse/J/physical gradient");
    Check(interface==0 && area==0 && minJ>0,"interface or global injectivity certificate");
    Check(wall<2e-8 && wallBound<=1e-9,"independent physical wall distance");
    Check(std::abs(volume-exact)<1e-11 && std::abs(volume-volume32)<1e-11 && std::abs(flux-2*volume)<1e-11,"independent volume/flux checks");
    Check(g.Signature()!=RANSTwoRegionHill(.04,.1,.02,.8,-.2,.2).Signature(),"geometry identity omitted Hw");
    g.RequireCompatible(g.Signature());
    Reject([&]{g.RequireCompatible(RANSTwoRegionHill(.04,.1,.02,.8,-.2,.2).Signature());},"incompatible geometry restore accepted");
    const RANSTwoRegionHill experimental(.04,.1,.015,.8,-1.,1.2);
    const auto crest=experimental.Locate({0,0,.05});
    Check(crest.inner && Norm(Sub(crest.chart,{0,0,.01}))<1e-12,"crest vertical .01m observation");
    const auto wake=experimental.Locate({.3,0,.01});
    Check(wake.inner && std::abs(wake.chart[0]-.3)<1e-12 && std::abs(wake.chart[2]-.01)<1e-12,"global x/Hh=7.5 observation");
    const auto plus=experimental.Locate(experimental.Inner(.03,.015).point,+1);
    const auto minus=experimental.Locate(experimental.Inner(.03,.015).point,-1);
    Check(!plus.inner && minus.inner,"interface observation side lost");
    Reject([&]{experimental.Locate({0,0,.02});},"solid observation accepted");
    const auto tied=g.PhysicalWallDistance({0,0,-.04},1e-9);
    Check(tied.footpoints.size()==2 && tied.footpoints[0][0]*tied.footpoints[1][0]<0,"multiple minimizing footpoints concealed");
    // A fluid-side tie between the flat floor and a hill footpoint. Neither
    // the scalar distance nor its consumers may select an averaged normal.
    const double fs=.02,dfs=-.2*std::acos(-1.),gbase=std::hypot(1.,dfs);
    const double radius=fs*gbase/(gbase-1),xTie=.05-fs*dfs/(gbase-1);
    const auto fluidTie=g.PhysicalWallDistance({xTie,0,radius},1e-9);
    Check(g.Locate({xTie,0,radius}).inner==false && fluidTie.footpoints.size()==2 && std::abs(fluidTie.distance-radius)<2e-9,
        "fluid-side hill/flat nearest-footpoint tie not retained");


    return J::Object({{"version",J::Quote("CosineNormalVerticalV1")},{"signature",J::Quote(g.Signature())},
        {"maximum_curvature",J::Number(g.CurvatureBound())},{"curvature_ratio",J::Number(.015*g.CurvatureBound())},
        {"minimum_sampled_outer_J",J::Number(minJ)},{"inverse_error",J::Number(inverse)},{"jacobian_fd_error",J::Number(jac)},
        {"physical_gradient_error",J::Number(gradient)},{"interface_coordinate_error",J::Number(interface)},
        {"opposite_area_vector_error",J::Number(area)},{"volume",J::Number(volume)},{"independent_graph_volume",J::Number(exact)},
        {"GL32_GL64_volume_difference",J::Number(std::abs(volume-volume32))},{"divergence_flux_error",J::Number(std::abs(flux-2*volume))},
        {"wall_distance_error",J::Number(wall)},{"wall_distance_interval_bound",J::Number(wallBound)},
        {"mapped_production_enabled","false"}});
}

std::string ReferenceChecks()
{
    using K=RANSFiniteInflowReference::Knot;std::vector<K> nodes;
    auto u=[](double z){return 1-std::pow(1-z/.3,3);};auto du=[](double z){return 10*std::pow(1-z/.3,2);};
    for(double z:{0.,.07,.2,.3})nodes.push_back({z,u(z),du(z),.02,.03});
    const RANSFiniteInflowReference f(nodes,{-.2,0,0},{1,0,0},{0,0,1},1.);
    double value=0,derivative=0,work=0,neutral=0,trace=0;
    for(double z:{.001,.04,.1,.22,.3,.8})
    {
        const auto v=f.AtHeight(z);value=std::max(value,std::abs(v.U-(z<=.3?u(z):1)));
        derivative=std::max(derivative,std::abs(v.dU-(z<=.3?du(z):0)));
        Check(v.gamma>=0 && v.k>0 && v.epsilon>0,"reference positivity");
        trace=std::max(trace,std::abs(v.strain[0][0]+v.strain[1][1]+v.strain[2][2]));
        RANSTensor g{};g[0][2]=v.dU;RANSRobustVariant a;a.beta=.4;a.gamma=v.gamma;a.referenceStrain=v.strain;
        auto b=a;b.gamma=0;const auto c=EvaluateRANSRobustClosure(RANSModel::CLSHighRe,v.k,v.epsilon,g,.1,{},0,a);
        const auto unregularized=EvaluateRANSRobustClosure(RANSModel::CLSHighRe,v.k,v.epsilon,g,.1,{},0,b);
        Check(c.reynoldsStress==unregularized.reynoldsStress,"gamma correction nonzero at prescribed strain");
        g[0][2]+=.1;const auto pert=EvaluateRANSRobustClosure(RANSModel::CLSHighRe,v.k,v.epsilon,g,.1,{},0,a);
        double p=0;for(int i=0;i<3;++i)for(int j=0;j<3;++j)p-=pert.reynoldsStress[i][j]*g[i][j];
        work=std::max(work,std::abs(p-pert.production));
        Check(std::abs(pert.epsilonProduction-1.44*v.epsilon/v.k*pert.production)<1e-12,"epsilon production stress-work");
        Check(f.AtHeight(z).gamma==v.gamma && f.AtHeight(z).strain==v.strain,"evolving fields modified reference");
    }
    const auto free=f.AtHeight(.8);Check(free.dU==0 && free.strain==RANSTensor{},"nonzero freestream reference shear");
    RANSTensor zero{};RANSRobustVariant v;v.beta=.4;v.gamma=free.gamma;v.referenceStrain=free.strain;
    const auto c=EvaluateRANSRobustClosure(RANSModel::CLSHighRe,free.k,free.epsilon,zero,.8,{},0,v);
    Check(c.reynoldsStress[0][2]==0 && c.production==0,"spurious free-stream residual shear");
    const auto cls=RANSNeutralReference::CLS();
    for(double z:{.01,.04,.1,1.})
    {
        const auto exact=cls.AtHeight(z,.074,.015,1.);
        // At an exact represented neutral knot the finite provider uses the
        // same Original viscosity, without imposing a separate log gamma law.
        RANSFiniteInflowReference n({{0,0,1,.01,.1},{z,exact.velocity,exact.velocityGradient,exact.k,exact.epsilon},{2*z,exact.velocity,0,exact.k,exact.epsilon}}, {},{1,0,0},{0,0,1},1.);
        const auto q=n.AtHeight(z);neutral=std::max(neutral,std::abs(q.gamma/(cls.Kappa()*.074*z)-1));
        Check(std::abs(q.strain[0][2]-exact.velocityGradient)<1e-12,"neutral strain recovery");
    }
    Check(value<1e-14 && derivative<1e-13 && neutral<1e-13 && work<1e-12 && trace==0,"finite reference exactness");
    double branchAction=0;
    {
        const auto f0=f.AtHeight(.1);RANSTensor g{};g[0][2]=f0.dU;
        RANSRobustVariant variant;variant.beta=.4;variant.gamma=f0.gamma;variant.referenceStrain=f0.strain;
        RANSTensor direction{};direction[0][2]=.2;direction[2][0]=.1;
        for(bool strain:{false,true})
        {
            auto eval=[&](double h)
            {auto v=g;for(int a=0;a<3;++a)for(int b=0;b<3;++b)v[a][b]+=h*direction[a][b];
             return EvaluateRANSRobustClosureLinearization(RANSModel::CLSHighRe,f0.k,f0.epsilon,v,.1,{},0,strain,variant);};
            const auto p=eval(1e-4),m=eval(-1e-4),p2=eval(5e-5),m2=eval(-5e-5);
            branchAction=std::max(branchAction,std::abs((p.production-m.production)/2e-4-(p2.production-m2.production)/1e-4));
        }
        Check(branchAction<1e-10,"fixed-branch directional action inconsistent");
        // The physical max(strain,rotation) selects strain for this positive
        // direction and rotation for its negative; do not central-difference
        // across the tie and call it a single classical derivative.
        for(double sign:{-1.,1.})
        {
            auto v=g;for(int a=0;a<3;++a)for(int b=0;b<3;++b)v[a][b]+=sign*1e-5*direction[a][b];
            const auto physical=EvaluateRANSRobustClosure(RANSModel::CLSHighRe,f0.k,f0.epsilon,v,.1,{},0,variant);
            const auto selected=EvaluateRANSRobustClosureLinearization(RANSModel::CLSHighRe,f0.k,f0.epsilon,v,.1,{},0,sign>0,variant);
            Check(physical.reynoldsStress==selected.reynoldsStress,"physical branch selection near strain/rotation tie");
        }
    }
    auto bad=nodes;bad.back().dU=1;Reject([&]{RANSFiniteInflowReference q(bad,{}, {1,0,0},{0,0,1},1);},"nonzero freestream derivative accepted");
    Reject([&]{f.AtHeight(-.001);},"negative incoming height accepted");
    return J::Object({{"version",J::Quote("FiniteInflowV1")},{"interpolation_error",J::Number(value)},
        {"derivative_error",J::Number(derivative)},{"neutral_viscosity_relative_error",J::Number(neutral)},
        {"stress_work_error",J::Number(work)},{"fixed_branch_action_difference",J::Number(branchAction)},{"trace_defect",J::Number(trace)},{"uniform_freestream_shear",J::Number(c.reynoldsStress[0][2])}});
}

std::string Preflight(const RANSFiniteInflowReference &profile,const std::string &queries)
{
    const auto cls=RANSNeutralReference::CLS();const auto wall=CalibratedAWFParameters(cls);
    const double nu=.04/10500,us=.051/.6867027875750893,primaryH=.015;
    const auto e=profile.AtHeight(primaryH);double a=0,b=primaryH*.999;
    // The legacy scalar AWF constructor needs an interior output coordinate;
    // J0 does not sample k there. Its K argument remains the endpoint trace.
    auto velocity=[&](double h){return us*us*AnalyticalWallFunction({nu,primaryH,h,primaryH/2},e.k,wall).Integrals().j0;};
    Check(velocity(a)>=e.U && velocity(b)<=e.U,"flat-only effective-roughness calibration not bracketed");
    for(int j=0;j<70;++j){double h=(a+b)/2;if(velocity(h)>e.U)a=h;else b=h;}
    const double h=(a+b)/2;
    std::vector<std::string> candidates;
    for(double H:{.01,.015,.02})
    {
        RANSTwoRegionHill geometry(.04,.1,H,.8,-.2,.2);const auto ref=profile.AtHeight(H);
        std::vector<std::string> scans;bool all=true;
        for(double curvature:{-geometry.CurvatureBound(),-geometry.CurvatureBound()/2,0.,geometry.CurvatureBound()/2,geometry.CurvatureBound()})for(double factor:{.5,1.,2.})
        {
            const double K=ref.k*factor;AWFProfileOptions options;options.formulation=AWFProfileFormulation::CLSCalibratedBudgetAWF;
            try
            {
                const EquilibriumAWFProfiles p({nu,H,h,H,H},K,cls,{curvature,0,.1},options,wall);
                const auto &v=p.Reference();Check(v.reason==AWFFallbackReason::Active && p.Dissipation(H)>0,"profile gate");
                std::vector<double> panels{0.,H,v.cutoff};
                for(double r:{h,p.Original().ViscosityInterceptHeight()})if(r>0 && r<H)panels.push_back(r);
                std::sort(panels.begin(),panels.end());panels.erase(std::unique(panels.begin(),panels.end()),panels.end());
                double independentE=0,independentP=0;
                for(std::size_t j=1;j<panels.size();++j)
                {
                    independentE+=Integrate([&](double r){return (1+curvature*r)*p.Dissipation(r);},panels[j-1],panels[j]);
                    independentP+=Integrate([&](double r){return (1+curvature*r)*p.ReferenceProduction(r);},panels[j-1],panels[j]);
                }
                const double independentError=std::abs(independentE-independentP)/independentP;
                Check(independentError<2e-11,"independent physical J-weighted cutoff budget");
                scans.push_back(J::Object({{"curvature",J::Number(curvature)},{"K",J::Number(K)},{"admissible","true"},
                    {"cutoff",J::Number(v.cutoff)},{"independent_J_budget_error",J::Number(independentError)},
                    {"counterfactual_UH",J::Number(v.matchingVelocity)},{"counterfactual_ustar",J::Number(v.frictionVelocity)},
                    {"counterfactual_z0",J::Number(H*std::exp(-cls.Kappa()*v.matchingVelocity/v.frictionVelocity))},{"budget_margin",J::Number(v.productionBudget/v.minimumBudget-1)},
                    {"relative_budget_error",J::Number(v.budgetResidual/v.productionBudget)},{"epsilon_H",J::Number(p.Dissipation(H))}}));
            }
            catch(const std::invalid_argument &error){all=false;scans.push_back(J::Object({{"curvature",J::Number(curvature)},{"K",J::Number(K)},{"admissible","false"},{"reason",J::Quote(error.what())}}));}
        }
        candidates.push_back(J::Object({{"Hw",J::Number(H)},{"h_eff",J::Number(h)},{"geometry_admissible","true"},
            {"curvature_ratio",J::Number(H*geometry.CurvatureBound())},{"K",J::Number(ref.k)},
            {"finite_UH",J::Number(ref.U)},{"finite_epsilon_H",J::Number(ref.epsilon)},
            {"K_range_factors",J::Numbers(std::vector<double>{.5,1.,2.})},
            {"profile_measurement_coverage",J::Quote(H<.01024086126?"slightly below first mean measurement: declared completion":"within permitted flat measured range")},
            {"all_declared_K_curvature_budget_cases_admissible",J::Boolean(all)},{"scan",J::List(scans)}}));
    }
    std::vector<std::string> tensors;double covariance=1e100;
    for(const auto &k:profile.Knots())if(k.z>=.01)
    {
        const auto p=profile.AtHeight(k.z);RANSTensor g{};g[0][2]=p.dU;
        RANSRobustVariant v;v.beta=.4;v.gamma=p.gamma;v.referenceStrain=p.strain;
        const auto c=EvaluateRANSRobustClosure(RANSModel::CLSHighRe,p.k,p.epsilon,g,k.z,{},0,v);
        const auto eig=J::Eigenvalues(c.reynoldsStress);covariance=std::min(covariance,eig[0]);
        tensors.push_back(J::Object({{"z",J::Number(k.z)},{"k",J::Number(p.k)},{"epsilon",J::Number(p.epsilon)},
            {"shear",J::Number(p.dU)},{"gamma",J::Number(p.gamma)},{"production",J::Number(c.production)},
            {"eigenvalues",J::Numbers(eig)},{"Rxx_Ryy_Rzz_Rxz",J::Numbers(std::vector<double>{c.reynoldsStress[0][0],c.reynoldsStress[1][1],c.reynoldsStress[2][2],c.reynoldsStress[0][2]})}}));
    }
    std::map<std::string,std::array<double,4>> errors;std::vector<std::string> queryResults;
    if(!queries.empty())
    {
        std::ifstream in(queries);double z,target,weight,bound;std::string quantity;
        while(in>>z>>quantity>>target>>weight>>bound)
        {
            const auto p=profile.AtHeight(z);RANSTensor g{};g[0][2]=p.dU;
            RANSRobustVariant v;v.beta=.4;v.gamma=p.gamma;v.referenceStrain=p.strain;
            const auto law=EvaluateRANSRobustClosure(RANSModel::CLSHighRe,p.k,p.epsilon,g,z,{},0,v);
            double predicted=quantity=="U"?p.U:(quantity=="minusRxz"?-law.reynoldsStress[0][2]:std::sqrt(law.reynoldsStress[0][0]));
            auto &e=errors[quantity];e[0]+=weight*std::pow(predicted-target,2);e[1]+=weight;e[2]+=weight*bound*bound;e[3]=std::max(e[3],std::abs(predicted-target));
            queryResults.push_back(J::Object({{"height",J::Number(z)},{"quantity",J::Quote(quantity)},{"target",J::Number(target)},
                {"completed_reference_prediction",J::Number(predicted)},{"weight",J::Number(weight)},{"extraction_bound",J::Number(bound)},
                {"inside_provisional_AWF_strip",J::Boolean(z<primaryH)}}));
        }
    }
    std::vector<std::pair<std::string,std::string>> fit;
    for(auto [name,e]:errors)fit.emplace_back(name,J::Object({{"weighted_RMSE",J::Number(std::sqrt(e[0]/e[1]))},
        {"retained_weight",J::Number(e[1])},{"weighted_extraction_envelope",J::Number(std::sqrt(e[2]/e[1]))},{"maximum_error",J::Number(e[3])}}));
    return J::Object({{"status",J::Quote("completed")},{"governing_solves","0"},
        {"roughness_calibration",J::Quote("provisional h_eff from permitted U_ref(.015), measured u_star ratio, completed endpoint K, calibrated AWF J0; no hill fit")},
        {"selected",J::Object({{"Hw",J::Number(primaryH)},{"h_eff",J::Number(h)},{"calibration_velocity_error",J::Number(velocity(h)-e.U)}})},
        {"flat_reference_fit_and_tensor_errors",J::Object(fit)},{"permitted_queries",J::List(queryResults)},
        {"tensor_comparison_scope",J::Quote("selected outer law on fixed completion; below Hw is a labelled counterfactual, actual pilot uses AWF there; no Fig3d or hill data fitted")},
        {"candidates",J::List(candidates)},{"finite_profile_covariance_minimum",J::Number(covariance)},{"profile_tensors",J::List(tensors)}});
}
}

int main(int argc,char **argv)
{
    std::string output="cao-h1-kernels.json",file,queries;
    for(int i=1;i+1<argc;i+=2){std::string k=argv[i];if(k=="--output")output=argv[i+1];else if(k=="--profile")file=argv[i+1];else if(k=="--queries")queries=argv[i+1];else return 2;}
    try
    {
        const auto reference=ReferenceChecks(),geometry=GeometryChecks();std::string preflight="null";
        if(!file.empty())
        {
            std::ifstream in(file);std::string text((std::istreambuf_iterator<char>(in)),{});
            const RANSFiniteInflowReference profile(RANSFiniteInflowReference::Parse(text),{-.2,0,0},{1,0,0},{0,0,1},1.);
            preflight=Preflight(profile,queries);
        }
        std::ofstream out(output);out<<J::Object({{"status",J::Quote("passed")},{"checks",J::Number(checks)},
            {"governing_solves","0"},{"physical_commits","0"},{"reference",reference},{"geometry",geometry},{"preflight",preflight}})<<'\n';
        std::cout<<"H1 native checks: "<<checks<<" passed\n";return 0;
    }
    catch(const std::exception&e){std::ofstream out(output);out<<J::Object({{"status",J::Quote("failed")},{"reason",J::Quote(e.what())},{"checks",J::Number(checks)}})<<'\n';std::cerr<<e.what()<<'\n';return 1;}
}
