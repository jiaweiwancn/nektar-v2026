// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeCoupling.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <SpatialDomains/MeshGraphIO.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <utility>

namespace
{
namespace fs=std::filesystem;
namespace J=Nektar::RANSWallSolveJson;
using namespace Nektar;
using namespace Nektar::WJNative;

void Require(bool condition,const std::string &message)
{
    if(!condition)throw std::runtime_error(message);
}

double Clock()
{
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

void VectorJSON(std::ostream &out,const Vec &v)
{
    out<<'[';
    for(Eigen::Index i=0;i<v.size();++i)out<<(i?",":"")<<v[i];
    out<<']';
}

std::unique_ptr<Coupling> MakeCoupling(const LibUtilities::SessionReaderSharedPtr &session,
                                      const SpatialDomains::MeshGraphSharedPtr &graph,
                                      InnerKTransport transport)
{
    std::string table;
    session->LoadSolverInfo("WJReferenceTable",table,"");
    auto reference=std::make_shared<Reference>(table,Reference::QualifiedTableSha256,transport);
    const auto &p=reference->Parameters();
    Geometry geometry;geometry.d0=p.d0;geometry.h=p.roughness;geometry.H=p.H;
    CoupledOptions options;
    options.sourceIdentity=WJ_INNER_K_FIXTURE_SOURCE_ID;
    options.backend=ConstitutiveBackend::FastFull3DV1;
    session->LoadParameter("TimeStep",options.dt,.05);
    session->LoadParameter("WJLinearTolerance",options.linearTolerance,1e-12);
    session->LoadParameter("WJCoupledTolerance",options.nonlinearTolerance,1e-10);
    session->LoadParameter("WJMaximumIterations",options.maximumIterations,60);
    session->LoadParameter("WJQuadratureX",options.quadratureX,8);
    session->LoadParameter("WJQuadratureR",options.quadratureR,12);
    session->LoadSolverInfo("WJInitialSourceIdentity",options.initialSourceIdentity,"");
    session->LoadSolverInfo("WJFixtureInitialReferenceIdentity",options.initialReferenceIdentity,"");
    session->LoadSolverInfo("WJSteadyVelocityPC",options.steadyVelocityPC,"FullGradientV1");
    std::vector<MultiRegions::ContFieldSharedPtr> fields(5);
    for(int i:{0,1,2,4})
        fields[i]=MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(session,graph,session->GetVariable(i),true,false);
    std::vector<unsigned int> outer;
    for(int i=0;i<fields[2]->GetExpSize();++i)
    {
        const auto exp=fields[2]->GetExp(i);
        double lo=1e100,hi=-1e100;
        for(int j=0;j<exp->GetGeom()->GetNumVerts();++j)
        {
            double x,z,d;
            exp->GetGeom()->GetVertex(j)->GetCoords(x,z,d);
            lo=std::min(lo,z-p.d0);hi=std::max(hi,z-p.d0);
        }
        Require(!(lo<p.H-1e-12&&hi>p.H+1e-12),"H must be a complete mesh edge");
        if(lo>=p.H-1e-12)outer.push_back(i);
    }
    Require(!outer.empty(),"empty outer epsilon mesh");
    auto subset=MemoryManager<MultiRegions::ExpList>::AllocateSharedPtr(*fields[2],outer,true);
    fields[3]=MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(session,*subset);
    return std::make_unique<Coupling>(fields,geometry,reference,options);
}

bool RejectHistory(Coupling &coupling,const fs::path &input)
{
    try{(void)coupling.ReadHistory(input);}
    catch(const std::exception &e)
    {return std::string(e.what()).find("checkpoint scientific identity differs")!=std::string::npos;}
    return false;
}

void ReferenceActions(std::ostream &out,const Reference &ref)
{
    const auto &p=ref.Parameters();
    bool comma=false;
    out<<'[';
    for(double K:{p.KMinimum,.947*p.KIn,p.KIn,1.063*p.KIn,p.KMaximum})
    {
        const double c=ref.Cutoff(K), deriv=ref.CutoffDerivative(K), inv=ref.KAtCutoff(c);
        for(double r:{0.,.7*c,1.3*c,.4*p.H,p.H})
        {
            const auto a=ref.Auxiliary(K,r), E=ref.DependentEpsilon(K,r);
            if(comma)out<<',';
            out<<"{\"K\":"<<K<<",\"r\":"<<r<<",\"c\":"<<c<<",\"c_K\":"<<deriv<<",\"inverse\":"<<inv
               <<",\"aux\":"<<a.value<<",\"aux_K\":"<<a.dK<<",\"aux_r\":"<<a.dr
               <<",\"E\":"<<E.value<<",\"E_K\":"<<E.dK<<",\"E_r\":"<<E.dr
               <<",\"D_inner\":"<<ref.KDiffusivity(K,E.value,a.value,false)
               <<",\"D_outer\":"<<ref.KDiffusivity(K,E.value,a.value,true)<<'}';
            comma=true;
        }
    }
    out<<']';
}
}

namespace Nektar::WJNative
{
struct InnerKFixtureAccess
{
    static CoupledHistory State(Coupling &c,bool moving)
    {
        auto state=c.Initialise().current;
        const auto &ref=c.ReferenceModel();const auto &p=ref.Parameters();
        double lower=0,upper=p.topRadius;
        for(const auto &e:c.m_flow.Elements())
        {
            if(e.rl<p.H-1e-12)lower=std::max(lower,e.rl);
            if(e.rr>p.H+1e-12)upper=std::min(upper,e.rr);
        }
        Require(lower<p.H && upper>p.H,"interface flux fixture needs both native neighbors");
        const double E=ref.Amplitude(p.KIn)/p.H;
        const double aux=ref.Aux(p.KIn,p.H);
        const double Dminus=ref.KDiffusivity(p.KIn,E,aux,false), Dplus=ref.KDiffusivity(p.KIn,E,aux,true);
        const double Q=.02*ref.Amplitude(p.KIn);
        state.k=c.m_k.Project([&](double x,double z){
            const double f=std::max(0.,1-std::abs(x)/2.5);
            const double delta=z-p.d0-p.H;
            if(moving)return p.KIn*(1+.04*f*(x-.30)/2.8);
            const double width=delta<0?p.H-lower:upper-p.H;
            const double t=std::abs(delta)/width;
            return p.KIn+(t<1?f*Q/(delta<0?Dminus:Dplus)*delta*(1-t)*(1-t):0.);
        });
        return c.FreshHistory(state);
    }

    static void Actions(std::ostream &out,Coupling &c,const CoupledHistory &h)
    {
        CoupledTiming timing;
        const auto data=c.Evaluate(h.current,0.,true,timing);
        const auto adv=c.Advection(h.current,data.kRules);
        const auto loads=c.Forcing(0.,timing,true);
        const auto current=c.CurrentActions(h.current,data);
        const auto diffusion=c.m_k.DiffusionAction(data.kRules,data.Dk,h.current.k);
        const Vec source=-c.KVolumeLoad(h.current,data);
        out<<"{\"k_diffusion\":";VectorJSON(out,diffusion);
        out<<",\"k_source\":";VectorJSON(out,source);
        out<<",\"k_transport\":";VectorJSON(out,adv.k);
        out<<",\"k_natural\":";VectorJSON(out,-loads.k);
        out<<",\"k_current\":";VectorJSON(out,current.k+adv.k-loads.k);
        out<<",\"epsilon_current\":";VectorJSON(out,current.epsilon+adv.epsilon-loads.epsilon);
        out<<",\"dependent_trace_error\":"<<data.traceError<<",\"cutoff_cells\":[";
        bool comma=false;
        for(std::size_t i=0;i<c.m_flow.Elements().size();++i)
        {
            const auto &e=c.m_flow.Elements()[i];
            const auto breaks=c.CurrentCutoffXBreaks(e,h.current);
            if(breaks.empty())continue;
            if(comma)out<<',';
            out<<"{\"owner\":"<<e.geometry<<",\"x_cuts\":"<<J::Numbers(breaks)
               <<",\"points\":"<<data.kRules[i]->w.size()<<'}';
            comma=true;
        }
        out<<"]}";
    }

    static void Flux(std::ostream &out,Coupling &c,const CoupledHistory &h)
    {
        const auto &ref=c.ReferenceModel();const auto &p=ref.Parameters();
        const double x=.45;
        out<<'[';
        bool comma=false;
        for(const auto &e:c.m_flow.Elements())
        {
            if(!(e.xl<x && x<e.xr))continue;
            const bool below=std::abs(e.rr-p.H)<1e-12,above=std::abs(e.rl-p.H)<1e-12;
            if(!below&&!above)continue;
            const auto rule=c.m_rules.HorizontalTrace(e,Vec::Constant(1,x),below?3:2);
            const double K=c.InterfaceK(e,x,h.current);
            const double k=c.m_k.Evaluate(rule,h.current.k)[0];
            const double E=above?c.m_e.Evaluate(rule,h.current.epsilon)[0]:ref.Ew(K,p.H);
            const double D=ref.KDiffusivity(k,E,below?ref.Aux(K,p.H):0.,above);
            const double grad=c.m_k.Evaluate(rule,h.current.k,1)[0];
            if(comma)out<<',';
            out<<"{\"owner\":"<<e.geometry<<",\"side\":"<<J::Quote(below?"lower":"upper")<<",\"x\":"<<x
               <<",\"K\":"<<K<<",\"E\":"<<E<<",\"D\":"<<D<<",\"gradient\":"<<grad<<",\"Q\":"<<D*grad<<'}';
            comma=true;
        }
        out<<']';
    }
};
}

namespace
{
void Verify(const LibUtilities::SessionReaderSharedPtr &session,const SpatialDomains::MeshGraphSharedPtr &graph,
            const fs::path &output)
{
    auto fixture=MakeCoupling(session,graph,InnerKTransport::AuxiliaryHistoricalFixtureV1);
    auto v2=MakeCoupling(session,graph,InnerKTransport::ScalarV2);
    const auto flux=InnerKFixtureAccess::State(*fixture,false);
    const auto moving=InnerKFixtureAccess::State(*fixture,true);
    const auto baseline=InnerKFixtureAccess::State(*v2,false);
    fixture->WriteHistory(output/"fixture.bin",flux);v2->WriteHistory(output/"v2.bin",baseline);
    const bool rejectV2=RejectHistory(*fixture,output/"v2.bin"),rejectFixture=RejectHistory(*v2,output/"fixture.bin");
    Require(rejectV2&&rejectFixture,"cross-transport history read accepted");
    std::ofstream out(output/"actions.json");out<<std::setprecision(17)
        <<"{\"schema\":\"WJInnerKFixtureNativeActionsV1\",\"cross_policy_histories_rejected\":true,\"TASK3_READY\":false"
        <<",\"kernel_source\":"<<J::Quote(WJ_H27_SOURCE_ID)<<",\"diagnostic_source\":"<<J::Quote(WJ_INNER_K_FIXTURE_SOURCE_ID)
        <<",\"reference_actions\":";
    ReferenceActions(out,fixture->ReferenceModel());
    out<<",\"H_flux\":";InnerKFixtureAccess::Flux(out,*fixture,flux);
    for(auto [name,history]:{std::make_pair("flux",&flux),std::make_pair("moving",&moving)})
    {
        out<<",\""<<name<<"\":{\"state\":";
        fixture->Report(out,*history,fixture->AuditSteadyState(history->current),true);
        out<<",\"actions\":";InnerKFixtureAccess::Actions(out,*fixture,*history);out<<'}';
    }
    out<<"}\n";Require(bool(out),"action output write failed");
}

void Steady(const LibUtilities::SessionReaderSharedPtr &session,const SpatialDomains::MeshGraphSharedPtr &graph,
            InnerKTransport transport,const fs::path &output)
{
    const double start=Clock();auto coupling=MakeCoupling(session,graph,transport);
    const double construction=Clock()-start;
    std::string initialFields;session->LoadSolverInfo("WJInitialFields",initialFields,"");
    const double init=Clock();
    auto history=initialFields.empty()?coupling->Initialise():coupling->ReadInitialFields(initialFields,true);
    const double initialization=Clock()-init;
    SteadyKrylovOptions policy;
    policy.acceptanceNorm=coupling->Options().nonlinearTolerance;
    session->LoadParameter("WJSteadySeconds",policy.maximumSeconds,0.);
    session->LoadParameter("WJSteadyNewtonIterations",policy.maximumNewtonIterations,50);
    session->LoadParameter("WJSteadyKrylovDimension",policy.restartDimension,80);
    session->LoadParameter("WJSteadyLinearCycles",policy.maximumLinearCycles,8);
    session->LoadParameter("WJSteadyFDStep",policy.fdRelativeStep,1e-6);
    Require(policy.maximumSeconds>0,"explicit steady solve wall cap required");
    CoupledStepInfo report;
    auto result=coupling->SolveSteady(history,policy,report);
    if(!result.converged)
    {
        std::ofstream fail(output/"failure.json");
        fail<<std::setprecision(17)<<"{\"status\":"<<J::Quote(result.status)<<",\"reason\":"<<J::Quote(result.reason)
            <<",\"residual\":"<<result.trueNorm<<",\"seconds\":"<<result.statistics.totalSeconds<<"}\n";
        throw std::runtime_error(result.reason);
    }
    Require(history.step==0 && history.time==0.,"steady fixture advanced time");
    coupling->WriteHistory(output/"fixture-steady.bin",history);
    std::ofstream out(output/"steady.json");out<<std::setprecision(17)
        <<"{\"schema\":\"WJH27SteadyCandidateV1\",\"normal_main\":false,\"diagnostic_only\":true"
        <<",\"qualified_reference\":false,\"TASK3_READY\":false,\"physical_time_advanced\":0,\"status\":\"CONVERGED\""
        <<",\"kernel_source\":"<<J::Quote(WJ_H27_SOURCE_ID)<<",\"diagnostic_source\":"<<J::Quote(WJ_INNER_K_FIXTURE_SOURCE_ID)
        <<",\"construction_seconds\":"<<construction<<",\"initialization_seconds\":"<<initialization
        <<",\"steady_seconds\":"<<result.statistics.totalSeconds<<",\"newton_steps\":"<<result.statistics.acceptedNewtonSteps
        <<",\"krylov_iterations\":"<<result.statistics.linearIterations<<",\"jacobian_actions\":"<<result.statistics.jacobianActions
        <<",\"true_residual\":"<<result.trueNorm<<",\"state\":";
    coupling->Report(out,history,report,true);out<<"}\n";Require(bool(out),"steady field write failed");
}
}

int main(int argc,char **argv)
{
    try
    {
        Require(argc==4,"usage: WJH27InnerKFixture actions|steady-v2|steady-fixture case.xml new-output-directory");
        const std::string mode=argv[1];
        Require(mode=="actions"||mode=="steady-v2"||mode=="steady-fixture","unknown diagnostic mode");
        const fs::path output=argv[3];
        Require(fs::create_directory(output),"diagnostic output directory exists");
        std::string xml=argv[2],opt=xml+".opt",program="WJH27InnerKFixture",option="--use-opt-file";
        char *args[]{program.data(),option.data(),opt.data(),xml.data(),nullptr};
        auto session=LibUtilities::SessionReader::CreateInstance(4,args);
        auto graph=SpatialDomains::MeshGraphIO::Read(session);
        std::string physical,solve,restart;
        int steps=0;
        session->LoadSolverInfo("WJNativeMode",physical,"");session->LoadSolverInfo("WJNativeSolve",solve,"");
        session->LoadSolverInfo("WJRestart",restart,"");session->LoadParameter("NumSteps",steps,0);
        Require(session->GetComm()->GetSize()==1 && physical=="Physical" && solve=="SteadyReference" && steps==0 && restart.empty(),
                "diagnostic requires one rank and a fresh unforced physical steady session");
        Require(session->GetVariables()==std::vector<std::string>({"u","v","k","epsilon","p"}),"wrong diagnostic fields");
        if(mode=="actions")Verify(session,graph,output);
        else Steady(session,graph,mode=="steady-v2"?InnerKTransport::ScalarV2:InnerKTransport::AuxiliaryHistoricalFixtureV1,output);
        session->Finalise();
        std::cout<<"INNER_K_FIXTURE_COMPLETED "<<mode<<'\n';
    }
    catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 2;}
}
