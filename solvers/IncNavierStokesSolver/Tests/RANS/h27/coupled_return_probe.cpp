// SPDX-License-Identifier: MIT
// Same public Coupling implementation as the registered application; a probe.
#include <IncNavierStokesSolver/EquationSystems/WJNativeCoupling.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <SpatialDomains/MeshGraphIO.h>
#include <boost/property_tree/json_parser.hpp>
#include <cstring>
#include <fstream>
#include <iostream>

using namespace Nektar;
using namespace Nektar::WJNative;
namespace J=RANSWallSolveJson;
namespace pt=boost::property_tree;
namespace fs=std::filesystem;
namespace
{
void Check(bool condition,const std::string&message)
{if(!condition)throw std::runtime_error(message);}
bool SameHistory(const CoupledHistory&a,const CoupledHistory&b)
{
    if(a.step!=b.step||a.time!=b.time||a.dt!=b.dt)return false;
    const std::array<const CoupledState*,4>left{&a.current,&a.previous,&a.advection,&a.previousAdvection},
        right{&b.current,&b.previous,&b.advection,&b.previousAdvection};
    for(int i=0;i<4;++i)for(const auto &pair:{std::make_pair(&left[i]->flow,&right[i]->flow),
        std::make_pair(&left[i]->k,&right[i]->k),std::make_pair(&left[i]->epsilon,&right[i]->epsilon)})
        if(pair.first->size()!=pair.second->size()||std::memcmp(pair.first->data(),pair.second->data(),pair.first->size()*sizeof(double)))return false;
    return true;
}
std::array<double, 6> Residuals(const CoupledResidual &r)
{
    return {r.momentum, r.continuity, r.gauge, r.k, r.epsilon, r.trace};
}
}
int main(int argc,char **argv)
{
    try
    {
        Check(argc==4,"usage: coupled_return_probe case.xml inputs.json new-output-directory");
        const fs::path output=argv[3];Check(!fs::exists(output)&&fs::create_directory(output),"immutable output directory exists");
        pt::ptree inputs;pt::read_json(argv[2],inputs);Check(inputs.get<std::string>("schema")=="WJH27CoupledReturnInputsV1","wrong input schema");
        std::string xml=argv[1],opt=xml+".opt",program="h27-coupled-return-probe",option="--use-opt-file";
        char *args[]{program.data(),option.data(),opt.data(),xml.data(),nullptr};
        auto session=LibUtilities::SessionReader::CreateInstance(4,args);auto graph=SpatialDomains::MeshGraphIO::Read(session);
        Check(session->GetComm()->GetSize()==1,"coupled return requires one rank");
        Check(session->GetVariables()==std::vector<std::string>({"u","v","k","epsilon","p"}),"wrong variable order");
        std::string table;session->LoadSolverInfo("WJReferenceTable",table,"");auto reference=std::make_shared<Reference>(table);
        Geometry geometry;geometry.d0=reference->Parameters().d0;geometry.H=reference->Parameters().H;geometry.h=reference->Parameters().roughness;
        CoupledOptions options;options.sourceIdentity=WJ_H27_SOURCE_ID;options.manufactured=true;options.checkCurrentActions=true;options.recordNonlinearTrace=true;
        session->LoadParameter("TimeStep",options.dt,.02);session->LoadParameter("WJLinearTolerance",options.linearTolerance,3e-13);
        session->LoadParameter("WJCoupledTolerance",options.nonlinearTolerance,3e-12);session->LoadParameter("WJMaximumIterations",options.maximumIterations,60);
        session->LoadParameter("WJQuadratureX",options.quadratureX,8);session->LoadParameter("WJQuadratureR",options.quadratureR,12);
        session->LoadParameter("WJMMSBodyAmplitude",options.mmsBodyAmplitude,1.);session->LoadSolverInfo("WJMMSVariant",options.mmsVariant,"Historical");
        options.backend=ConstitutiveBackend::FastFull3DV1;options.advancement="CoupledJFNKV1";options.stepPreconditioner="CurrentDataDriftV1";
        session->LoadParameter("WJStepPCRelativeDrift",options.stepPCRelativeDrift,.02);session->LoadParameter("WJStepPCMaximumUses",options.stepPCMaximumUses,8);
        options.stepNewton.acceptanceNorm=options.nonlinearTolerance;
        session->LoadParameter("WJStepNewtonSeconds",options.stepNewton.maximumSeconds,30.);
        session->LoadParameter("WJStepNewtonIterations",options.stepNewton.maximumNewtonIterations,20);
        session->LoadParameter("WJStepKrylovDimension",options.stepNewton.restartDimension,40);
        session->LoadParameter("WJStepLinearCycles",options.stepNewton.maximumLinearCycles,8);
        session->LoadParameter("WJStepFDStep",options.stepNewton.fdRelativeStep,1e-6);
        auto make=[&](CoupledOptions policy)
        {
            std::vector<MultiRegions::ContFieldSharedPtr>fields(5);
            for(int i:{0,1,2,4})fields[i]=MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(session,graph,session->GetVariable(i),true,false);
            std::vector<unsigned int>outer;
            for(int i=0;i<fields[2]->GetExpSize();++i)
            {
                const auto exp=fields[2]->GetExp(i);double low=1e100,high=-1e100;
                for(int j=0;j<exp->GetGeom()->GetNumVerts();++j){double x,z,d;exp->GetGeom()->GetVertex(j)->GetCoords(x,z,d);low=std::min(low,z-geometry.d0);high=std::max(high,z-geometry.d0);}
                Check(!(low<geometry.H-1e-12&&high>geometry.H+1e-12),"H is not a complete edge layer");if(low>=geometry.H-1e-12)outer.push_back(i);
            }
            auto subset=MemoryManager<MultiRegions::ExpList>::AllocateSharedPtr(*fields[2],outer,true);
            fields[3]=MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(session,*subset);
            return std::make_unique<Coupling>(fields,geometry,reference,policy);
        };
        auto warm=make(options);
        const auto aPath=inputs.get<std::string>("A.path");
        const auto bPath=inputs.get<std::string>("B.path");
        const auto A=warm->ReadHistory(aPath), B=warm->ReadHistory(bPath);
        Check(A.step==10&&A.time==.2&&B.step==18&&B.time==.36,"wrong actual A/B histories");
        int sequence=0;
        auto replay=[&](const char *name,Coupling &coupling,const CoupledHistory &input,
                        int instance,bool reset,bool retry,bool expectFailure)
        {
            const fs::path stem=output/name;auto history=input;coupling.WriteHistory(stem.string()+".input.bin",history);
            if(reset)coupling.ResetPreconditioners();
            // At the exact supplied coefficients and clock, exercise the
            // complete zero-time-mass action with current advection and the
            // analytic time derivative removed from the manufactured load.
            const auto audit=coupling.AuditSteadyState(input.current,input.time);
            CoupledStepInfo info;bool accepted=false;std::string failure;
            try{info=coupling.Advance(history,retry);accepted=true;}
            catch(const std::exception &error){failure=error.what();info=coupling.LastAttemptInfo();}
            Check(accepted!=expectFailure,"unexpected acceptance/rejection in "+std::string(name)+"; "+failure);
            const bool unchanged=SameHistory(history,input);
            if(expectFailure)Check(unchanged&&failure.find("NEWTON_NOT_CONVERGED")!=std::string::npos&&info.retries==0,"wrong failure or accepted-history mutation");
            else Check(history.step==input.step+1&&info.residual.Maximum()<=options.nonlinearTolerance,"current step failed acceptance");
            coupling.WriteHistory(stem.string()+".final.bin",history);
            std::ofstream file(stem.string()+".result.json");file<<std::setprecision(17)
                <<"{\"schema\":\"WJH27CoupledReturnCaseV1\",\"normal_main\":false,\"actual_native_Coupling\":true,\"case\":"<<J::Quote(name)
                <<",\"sequence\":"<<sequence++<<",\"instance\":"<<instance<<",\"explicit_reset\":"<<(reset?"true":"false")
                <<",\"injected_retry\":"<<(retry?"true":"false")<<",\"accepted\":"<<(accepted?"true":"false")
                <<",\"history_unchanged\":"<<(unchanged?"true":"false")<<",\"failure_reason\":"<<J::Quote(failure)
                <<",\"exact_input_current_audit\":"<<J::Object({
                    {"time",J::Number(input.time)},{"step",J::Number(input.step)},
                    {"residuals",J::Numbers(Residuals(audit.residual))},{"seconds",J::Number(audit.seconds)},
                    {"minimum_k",J::Number(audit.minimumK)},{"minimum_epsilon",J::Number(audit.minimumEpsilon)},
                    {"minimum_covariance_over_k",J::Number(audit.minimumCovarianceOverK)},
                    {"minimum_principal_over_mu",J::Number(audit.minimumPrincipalOverMu)}})
                <<",\"state\":";
            coupling.Report(file,history,info,true);file<<"}\n";Check(bool(file),"result write failed");
        };
        replay("warm-a-first",*warm,A,0,false,false,false);
        replay("warm-b",*warm,B,0,false,false,false);
        std::string staleReason;
        try{(void)warm->ReadHistory(inputs.get<std::string>("stale_A.path"));}
        catch(const std::exception &error){staleReason=error.what();}
        Check(staleReason.find("restored explicit history differs")!=std::string::npos,
              "stale explicit history was not rejected by native reconstruction");
        {
            std::ofstream file(output/"stale-rejection.json");
            file<<J::Object({{"schema",J::Quote("WJH27StaleExplicitRejectionV1")},
                             {"source_path",J::Quote(inputs.get<std::string>("stale_A.path"))},
                             {"reason",J::Quote(staleReason)}})<<'\n';
            Check(bool(file),"stale-history result write failed");
        }
        const auto returned=warm->ReadHistory(aPath);
        Check(SameHistory(returned,A),"A did not reconstruct to identical bytes");
        replay("warm-a-return",*warm,returned,0,false,false,false);
        replay("retry-a",*warm,returned,0,false,true,false);
        auto cold=make(options);
        replay("cold-a",*cold,cold->ReadHistory(aPath),1,true,false,false);
        replay("cold-b",*cold,cold->ReadHistory(bPath),1,true,false,false);
        auto limited=options;limited.stepNewton.maximumNewtonIterations=1;auto failure=make(limited);
        replay("failed-a",*failure,A,2,false,false,true);
        session->Finalise();std::cout<<"SOURCE19_COUPLED_RETURN_REPLAYS_RECORDED\n";
    }
    catch(const std::exception &error){std::cerr<<error.what()<<'\n';return 2;}
}
