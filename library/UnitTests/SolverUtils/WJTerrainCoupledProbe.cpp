// SPDX-License-Identifier: MIT
// A small native coupled verification driver, not an experimental inlet/outlet.
// Governing steps call the production patch and the native IMEX history code.
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallTimeIntegration.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionResolution.h>
#include <LibUtilities/Communication/CommSerial.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <SpatialDomains/MeshGraphIO.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
using namespace Nektar;
using namespace Nektar::SolverUtils;
using Vector=RANSWallAlgebra::Vector;
using Vec=RANSVector;
using Tensor=RANSTensor;
using RANSWallAlgebra::Require;
namespace J=Nektar::RANSWallSolveJson;
namespace
{
double Dot(const Vec &a,const Vec &b){double r=0;for(int i=0;i<3;++i)r+=a[i]*b[i];return r;}
double Volume(const RANSWallEvolutionPatch &p)
{double v=0;for(const auto &q:p.Quadrature(0))v+=q.weight;return v;}
std::string TensorJSON(const Tensor &g)
{std::vector<std::string> rows;for(const auto &r:g)rows.push_back(J::Numbers(r));return J::List(rows);}
std::string AuditJSON(const RANSWallEvolutionPatch::Diagnostic &d)
{std::vector<std::pair<std::string,std::string>> rows;for(const auto &[k,v]:d)rows.push_back({k,J::Numbers(v)});return J::Object(rows);}
void Save(const std::string &path,const std::string &value)
{std::ofstream f(path);f<<value<<'\n';f.close();Require(f.good(),"terrain packet write failed");}
Vec Velocity(const RANSWallEvolutionPatch &p,const Vector &s,const RANSWallGlobalPoint &q,Tensor &g)
{Vec u{};for(int a=0;a<p.Dimension();++a){const int aa=p.PhysicalAxis(a);u[aa]=q.Value(s,p.offset[a]);g[aa]=q.Gradient(s,p.offset[a]);}return u;}

std::string Interface(const RANSWallEvolutionPatch &p,const Vector &s,double time,const RANSWallEvolutionData &data)
{
    std::vector<std::string> rows;double common=0;
    for(const auto &f:p.Surface(0,2,p.Parameters().geometry.height,1,12))
    {
        const auto q=f.test.coordinate;const auto frame=p.Mapping().At(q);
        const auto vi=p.Point(0,q,-1),vo=p.Point(0,q,1),ep=p.Point(p.Epsilon(),q,1);
        const auto ki=vi.Gradient(s,p.offset[p.K()]),ko=vo.Gradient(s,p.offset[p.K()]);
        Tensor gi{},go{},strain{};const auto ui=Velocity(p,s,vi,gi),uo=Velocity(p,s,vo,go);
        for(int a=0;a<3;++a)for(int b=0;b<3;++b)strain[a][b]=gi[a][b]+gi[b][a];
        const double k=vi.Value(s,p.offset[p.K()]),e=ep.Value(s,p.offset[p.Epsilon()]);
        const auto profile=p.StateProfile(s,q,time,data);const double nt=profile->Original().EddyViscosity(frame.height);
        const auto inner=p.InnerStress(k,nt,gi);
        const auto outer=p.OuterClosure(k,e,go,q);
        const double pi=p.Point(p.InnerPressure(),q,-1).Value(s,p.offset[p.InnerPressure()]);
        const double po=p.Point(p.OuterPressure(),q,1).Value(s,p.offset[p.OuterPressure()]);
        const double mechanicalI=WJMechanicalPressure(pi,k,0),mechanicalO=WJMechanicalPressure(po,k,0);
        const auto stress=RANSWallInnerMeanStress(inner,p.Parameters().geometry.molecularNu,pi);
        const auto to=WJOuterTraction(outer.reynoldsStress,mechanicalO,frame.normal);
        Vec ti{},delta{},ku{},kother{};
        for(int a=0;a<3;++a)
        {
            ti[a]=Dot(stress[a],frame.normal);delta[a]=ti[a]-to[a];
            ku[a]=ui[a]*k-(p.Parameters().geometry.molecularNu+nt/p.SigmaK())*ki[a];
            kother[a]=uo[a]*k-outer.kDiffusion*ko[a];
        }
        const double jump=WJInterfacePressureJump(k,p.Parameters().geometry.molecularNu+nt,strain,outer.reynoldsStress,frame.normal);
        Require(std::abs(Dot(delta,frame.normal)-(jump-mechanicalI+mechanicalO))<1e-10,
            "normal pressure-jump/full-traction identity differs");
        Require(vi.ids==vo.ids,"shared interface velocity map differs");
        for(int j=0;j<vi.ids.size();++j)common=std::max(common,std::abs(vi.phi[j]-vo.phi[j]));
        rows.push_back(J::Object({{"chart",J::Numbers(q)},{"physical",J::Numbers(p.Physical(q))},
            {"inner_geometry",J::Number(vi.geometry)},{"outer_geometry",J::Number(vo.geometry)},
            {"surface_weight",J::Number(f.weight*std::hypot(f.area[0],f.area[1],f.area[2]))},
            {"normal",J::Numbers(frame.normal)},{"tangent",J::Numbers(frame.tangent1)},
            {"K",J::Number(profile->K())},{"k",J::Number(k)},{"epsilon",J::Number(e)},
            {"epsilon_target",J::Number(profile->Dissipation(frame.height))},
            {"cutoff",J::Number(profile->Reference().cutoff)},{"cutoff_gate",J::Number(profile->Reference().productionBudget/profile->Reference().minimumBudget)},
            {"metric_stretch",J::Number(frame.J)},{"traction_inner",J::Numbers(ti)},{"traction_outer",J::Numbers(to)},
            {"traction_defect",J::Numbers(delta)},{"mechanical_jump",J::Number(mechanicalI-mechanicalO)},
            {"normal_traction_jump_target",J::Number(jump)},{"mass_inner",J::Number(Dot(ui,frame.normal))},
            {"mass_outer",J::Number(-Dot(uo,frame.normal))},{"k_flux_inner",J::Number(Dot(ku,frame.normal))},
            {"k_flux_outer",J::Number(-Dot(kother,frame.normal))},{"G_inner",TensorJSON(gi)},{"G_outer",TensorJSON(go)},
            {"R_inner",TensorJSON(inner.reynoldsStress)},{"R_outer",TensorJSON(outer.reynoldsStress)},
            {"production_inner",J::Number(inner.production)},{"production_outer",J::Number(outer.production)},
            {"nu_over_nut_H",J::Number(p.Parameters().geometry.molecularNu/nt)}}));
    }
    return J::Object({{"samples",J::List(rows)},{"shared_value_test_row_defect",J::Number(common)},
        {"scope",J::Quote("shared trace tests give opposite common mass/advection/traction loads; one-sided physical traction and k-flux defects are separately retained, not imposed as extra conditions")}});
}

// New independent quadrature on the actual represented fields. Only three
// declared test modes are exported; their native full coefficient contractions
// are retained. This is not an all-mode or nonlinear-Jacobian certificate.
std::string Weak(const RANSWallEvolutionPatch &p,const Vector &s,const Vector &h,double alpha,double time,
    const RANSWallEvolutionData &data)
{
    const auto hash=RANSWallR3::Hash(s);const auto work=p.InitializationWork();
    const double H=p.Parameters().geometry.height;
    std::vector<std::vector<Vector>> tests(3);
    for(int mode=0;mode<3;++mode)for(int block=0;block<p.Gauge();++block)
        tests[mode].push_back(p.fields[block]->Project([&](const Vec &q){
            if(mode==0)return 1.;if(mode==1)return q[0]/5;
            return q[2]<H/2 || q[2]>2?0.:q[2]<H?(q[2]-H/2)/(H/2):(2-q[2])/(2-H);}));
    auto contract=[&](const Vector &r)
    {
        std::vector<std::string> modes;
        for(int mode=0;mode<3;++mode)
        {Vector values;for(int b=0;b<p.Gauge();++b){long double v=0;
         for(int j=0;j<p.size[b];++j)v+=(long double)tests[mode][b][j]*r[p.offset[b]+j];values.push_back(v);}
         modes.push_back(J::Numbers(values));}return J::List(modes);
    };
    const auto components=p.ResidualResolutionComponents(s,h,alpha,time,data,false);
    std::vector<std::pair<std::string,std::string>> terms;
    for(const auto &[name,values]:components.terms)terms.push_back({name,contract(values)});
    const auto inner=p.Residual(s,h,alpha,time,data,false,false,0),outer=p.Residual(s,h,alpha,time,data,false,false,1);
    std::vector<std::string> rows;
    auto append=[&](const RANSWallGlobalPoint &v,const RANSWallGlobalPoint &pressure,const RANSWallGlobalPoint &ep,
        double weight,bool in,const EquilibriumAWFProfiles *profile)
    {
        Tensor G{};const auto u=Velocity(p,s,v,G);const auto kgrad=v.Gradient(s,p.offset[p.K()]);
        const auto ge=in?Vec{}:ep.Gradient(s,p.offset[p.Epsilon()]);
        const int pb=in?p.InnerPressure():p.OuterPressure();
        const double k=v.Value(s,p.offset[p.K()]),e=in?profile->Dissipation(v.coordinate[2]):ep.Value(s,p.offset[p.Epsilon()]);
        std::vector<std::string> modes;
        for(int mode=0;mode<3;++mode)
        {const auto tg=v.Gradient(tests[mode][0]);const auto eg=in?Vec{}:ep.Gradient(tests[mode][p.Epsilon()]);
         modes.push_back(J::Object({{"v",J::Number(v.Value(tests[mode][0]))},{"grad_v",J::Numbers(tg)},
          {"p",J::Number(pressure.Value(tests[mode][pb]))},{"e",in?"0":J::Number(ep.Value(tests[mode][p.Epsilon()]))},{"grad_e",J::Numbers(eg)}}));}
        const auto metric=p.Mapping().At(v.coordinate);
        rows.push_back(J::Object({{"inner",J::Boolean(in)},{"chart",J::Numbers(v.coordinate)},
            {"physical",J::Numbers(p.Physical(v.coordinate))},{"geometry",J::Number(v.geometry)},
            {"element",J::Number(v.element)},{"weight",J::Number(weight)},{"U",J::Numbers(u)},{"G",TensorJSON(G)},
            {"k",J::Number(k)},{"epsilon",J::Number(e)},{"grad_k",J::Numbers(kgrad)},{"grad_epsilon",J::Numbers(ge)},
            {"pi",J::Number(pressure.Value(s,p.offset[pb]))},{"lambda",J::Number(s[p.offset[p.Gauge()]])},
            {"history",J::Numbers(Vector{v.Value(h,p.offset[0]),v.Value(h,p.offset[1]),v.Value(h,p.offset[p.K()]),in?0:ep.Value(h,p.offset[p.Epsilon()])})},
            {"K",in?J::Number(profile->K()):"null"},{"cutoff",in?J::Number(profile->Reference().cutoff):"null"},
            {"kappa_eff",in?J::Number(metric.kappa1+metric.kappa2):"null"},
            {"force",J::Numbers(in?p.CartesianForcing(s,v.coordinate,time,data,false):Vec{})},{"modes",J::List(modes)}}));
    };
    const auto rule=LibUtilities::PointsManager()[LibUtilities::PointsKey(24,LibUtilities::eGaussGaussLegendre)];
    for(const auto &el:p.innerVelocity->Elements())if(el.lower[2]==0)
        for(int j=0;j<24;++j)
        {
            Vec anchor{(el.lower[0]+el.upper[0])/2+rule->GetZ()[j]*(el.upper[0]-el.lower[0])/2,0,0};
            const auto selected=p.StateProfile(s,anchor,time,data);auto options=selected->Options();options.quadratureOrder=32;
            const auto metric=p.Mapping().At(anchor);
            auto profile=EquilibriumAWFProfiles::WJ(selected->Geometry(),selected->K(),p.Parameters().model,options,
                {metric.kappa1,metric.kappa2});
            const double wx=rule->GetW()[j]*(el.upper[0]-el.lower[0])/2*metric.baseArea;
            for(const auto &q:profile.Quadrature())
            {auto a=anchor;a[2]=q.r;append(p.Point(0,a),p.Point(p.InnerPressure(),a),{},wx*q.weight*q.jacobian,true,&profile);}
        }
    for(const auto &v:p.outerVelocity->Quadrature(24))
    {
        const auto a=v.coordinate;
        append(p.Point(0,a),p.Point(p.OuterPressure(),a),p.Point(p.Epsilon(),a),v.weight*p.Mapping().At(a).determinant,false,nullptr);
    }
    Require(hash==RANSWallR3::Hash(s) && work==p.InitializationWork(),"weak observation mutated state/Lift");
    return J::Object({{"state_sha256",J::Quote(hash)},{"alpha_mass",J::Number(alpha)},
        {"scope",J::Quote("three physical weak modes:1,s/5,interface_hat; raw regional loads retain constrained exterior reactions; no additional epsilon flux")},
        {"native_inner",contract(inner)},{"native_outer",contract(outer)},
        {"native_terms",J::Object(terms)},{"volume",J::Number(Volume(p))},{"samples",J::List(rows)},
        {"governing_solves","0"}});
}
}
int main(int argc,char **argv)
{
    try
    {
        Require(argc==3,"usage: WJTerrainCoupledProbe SESSION.xml OUTPUT.json");
        const auto started=std::chrono::steady_clock::now();
        auto world=LibUtilities::SessionReader::CreateInstance(2,argv);world->InitSession();
        char program[]="wj-terrain-coupled";char *args[]{program,nullptr};
        auto session=LibUtilities::SessionReader::CreateInstance(1,args,std::vector<std::string>{argv[1]},LibUtilities::CommSerial::create(1,args));
        session->InitSession();session->SetSolverInfo("GlobalSysSoln","DirectFull");
        const auto graph=SpatialDomains::MeshGraphIO::Read(session);RANSWallPatchParameters par;
        session->LoadParameter("RANSWallHeight",par.geometry.height,.375);session->LoadParameter("RANSWallRoughness",par.geometry.roughness,.09);
        session->LoadParameter("Kinvis",par.geometry.molecularNu,.051*.005/3.97);par.geometry.sampleHeight=par.geometry.tailStart=par.geometry.height;par.chi=0;
        session->LoadParameter("RANSWallTangentialOrder",par.tangentialOrder,16);session->LoadParameter("RANSWallNormalOrder",par.normalOrder,24);
        LoadWJNativeConfiguration(session,par);
        RANSWallEvolutionPatch p(session,graph,world->GetComm(),par,LoadRANSWallNormalMap(session,2,par.geometry.height));
        Require(p.WJTerrain() && p.Size()<=3000,"terrain prototype requires its declared small coupled route");
        double us,dt;int steps;session->LoadParameter("RANSWallFrictionVelocity",us,.051);
        session->LoadParameter("TimeStep",dt,1e-6);session->LoadParameter("NumSteps",steps,3);
        Require(steps>=0 && steps<=3 && dt>0,"terrain prototype step/time bound");
        const double H=par.geometry.height,K=par.model.Q()*us*us,E=us*us*us/(par.model.Kappa()*H);
        std::vector<RANSWallEvolutionFunction> functions;
        for(int b=0;b<p.Gauge();++b)functions.push_back([&,b](const Vec &q,double){
            const double t=q[2]/H;
            if(b<2){if(q[2]>=H)return 0.;const auto m=p.Mapping().At(q);
                const double speed=.1*2*t*(1-t)*(1-2*t);return speed*m.tangent1[p.PhysicalAxis(b)];}
            if(b==p.K())return q[2]>=H?K:K*(1+.1*(1-3*t*t+2*t*t*t));
            if(b==p.Epsilon())return E;return 0.;});
        auto state=p.Project(functions,0.);
        // Declare initial internal velocity projection exactly at the zero
        // streamfunction H trace. This is initialization only, not an evolving
        // interface condition or a step restoration.
        for(int b=0;b<2;++b)for(auto [id,value]:p.fields[b]->Boundary([](const Vec &){return 0.;},
            [&](int axis,int side,const RANSWallGlobalField::Element &el){return axis==2 &&
                ((side && el.upper[2]==H) || (!side && el.lower[2]==H));}))
            { (void)value;state[p.offset[b]+id]=0.; }
        RANSWallEvolutionData data(2);data.fixedBoundary=p.FixedReferenceBoundary(state);
        state=p.Lift(state,0,data);p.SetReference(state,data,0,0);
        const auto initial=state;const auto initialWork=p.InitializationWork();
        const auto signature=std::string("WJH1-VerificationV1:")+p.Signature()+":input="+RANSWallR3::HashFile(argv[1])+":dt="+J::Number(dt);
        const std::string base=std::string(argv[2])+".checkpoints";
        std::vector<std::string> events,callbacks,checkpoints;Vector lastHistory=state;double lastAlpha=0;
        auto startup=p.CompleteRestrictedAudit(state,state,0,0,data,true);
        Save(std::string(argv[2])+".startup-audit.json",AuditJSON(startup));
        Require(startup.at("wj_hard_admissible")[0]==1,"terrain startup hard model-domain audit failed");
        const auto startInterface=Interface(p,state,0,data);
        const auto startAssembly=p.WJTerrainAssemblyChecks(state,state,1/dt,0,data);
        Save(std::string(argv[2])+".startup-weak.json",Weak(p,state,state,0,0,data));
        RANSWallTimeIntegration integration(1,dt,
            [&](const Vector &h,double time,double lambda)
            {
                auto history=state;std::copy(h.begin(),h.end(),history.begin());int iterations=0;
                auto candidate=p.Solve(state,history,1/lambda,time,data,false,iterations);
                Require(p.PimpleRecords().back().find("\"bound_activations\":0")!=std::string::npos,
                    "terrain StopOnBound: material numerical correction stops before commit");
                const auto audit=p.CompleteRestrictedAudit(candidate,history,1/lambda,time,data,true);
                Save(std::string(argv[2])+".candidate.json",J::Object({{"state",J::Numbers(candidate)},
                    {"history",J::Numbers(history)},{"time",J::Number(time)},{"alpha",J::Number(1/lambda)},{"audit",AuditJSON(audit)}}));
                Require(audit.at("wj_hard_admissible")[0]==1,"terrain completed candidate hard model-domain audit failed");
                callbacks.push_back(J::Object({{"state",J::Numbers(candidate)},{"history",J::Numbers(history)},
                    {"state_sha256",J::Quote(RANSWallR3::Hash(candidate))},{"history_sha256",J::Quote(RANSWallR3::Hash(history))},
                    {"time",J::Number(time)},{"alpha",J::Number(1/lambda)},{"audit",AuditJSON(audit)}}));
                lastHistory=history;lastAlpha=1/lambda;state=std::move(candidate);
                return Vector(state.begin(),state.begin()+p.DifferentialSize());
            },[&](const Vector &h,double time){auto s=state;std::copy(h.begin(),h.end(),s.begin());return p.ExplicitBody(s,time,data);});
        integration.Initialise(Vector(state.begin(),state.begin()+p.DifferentialSize()),0);
        auto saveCheckpoint=[&](){RANSWallTimeIntegration::Checkpoint cp;cp.signature=signature;cp.integration=integration.Capture();
            cp.coupledState=state;cp.fixedReference=initial;cp.fixedDefect=p.Defect();cp.eta=0;
            checkpoints.push_back(J::Quote(RANSWallTimeIntegration::WriteCheckpoint(base,cp)));};
        saveCheckpoint();
        for(int j=0;j<steps;++j)
        {
            integration.Advance();saveCheckpoint();
            events.push_back(J::Object({{"step",J::Number(integration.Steps())},{"time",J::Number(integration.Time())},
                {"interface",Interface(p,state,integration.Time(),data)},{"physical_samples",AuditJSON(p.PhysicalSamples(state,integration.Time(),data))}}));
        }
        const auto weak=std::string(argv[2])+".final-weak.json";
        Save(weak,Weak(p,state,lastHistory,lastAlpha,integration.Time(),data));
        Save(argv[2],J::Object({{"schema",J::Quote("WJH1NativeCoupledV1")},{"signature",J::Quote(signature)},
            {"map",J::Quote(p.Mapping().Signature())},{"state_size",J::Number(p.Size())},{"field_sizes",J::Numbers(p.size)},
            {"differential_size",J::Number(p.DifferentialSize())},{"volume",J::Number(Volume(p))},
            {"steps",J::Number(integration.Steps())},{"time",J::Number(integration.Time())},{"dt",J::Number(dt)},
            {"initial",J::Numbers(initial)},{"state",J::Numbers(state)},{"initialization_work",J::Numbers(initialWork)},
            {"startup_audit",AuditJSON(startup)},{"startup_interface",startInterface},{"startup_assembly_checks",startAssembly},
            {"final_assembly_checks",p.WJTerrainAssemblyChecks(state,lastHistory,steps?lastAlpha:1/dt,integration.Time(),data)},
            {"callbacks",J::List(callbacks)},{"events",J::List(events)},{"pimple_attempts",J::List(p.PimpleRecords())},
            {"checkpoints",J::List(checkpoints)},{"final_weak",J::Quote(weak)},
            {"roughness_budget",J::Numbers(p.RoughnessBudget(state,integration.Time(),data))},
            {"energy_budget",J::Numbers(p.EnergyBudget(state,lastHistory,lastAlpha,integration.Time(),data))},
            {"seconds",J::Number(std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count())},
            {"scope",J::Quote("controlled streamfunction vortex with reference-Dirichlet exterior; production coupled patch/native IMEX, not calibrated experimental hill inflow or open boundaries")}}));
        session->Finalise();world->Finalise();return 0;
    }
    catch(const std::exception &error){std::cerr<<error.what()<<std::endl;return 1;}
}
