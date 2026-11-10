// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallTimeIntegration.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallPhysicsSignature.h>
#include <cstring>
#include <LibUtilities/Communication/CommSerial.h>
#include <SpatialDomains/MeshGraphIO.h>
#include <boost/property_tree/json_parser.hpp>
#include <chrono>
#include <fstream>
#include <iostream>
using namespace Nektar;
using namespace Nektar::SolverUtils;
namespace J=Nektar::RANSWallSolveJson;
using Tree=boost::property_tree::ptree;
using Vector=RANSWallAlgebra::Vector;
namespace
{
Vector Values(const Tree &p,const std::string &key)
{Vector v;for(const auto &x:p.get_child(key))v.push_back(x.second.get_value<double>());return v;}
void Require(bool ok,const std::string &message){if(!ok)throw std::runtime_error(message);}
}
int main(int argc,char **argv)
{
    try
    {
        Require(argc==4,"usage: RANSPreservationProbe SESSION.xml REQUEST.json RESPONSE.json");
        const auto start=std::chrono::steady_clock::now();
        Tree request;boost::property_tree::read_json(argv[2],request);
        auto world=LibUtilities::SessionReader::CreateInstance(2,argv);world->InitSession();
        Require(world->GetComm()->GetSize()==1,"preservation verification requires serial execution");
        char name[]="preservation-probe";char *args[]{name,nullptr};
        std::vector<std::string> files{argv[1]};
        auto session=LibUtilities::SessionReader::CreateInstance(1,args,files,LibUtilities::CommSerial::create(1,args));session->InitSession();
        session->SetSolverInfo("GlobalSysSoln","DirectFull");auto graph=SpatialDomains::MeshGraphIO::Read(session);
        RANSWallPatchParameters p;
        session->LoadParameter("RANSWallChi",p.chi,1.);session->LoadParameter("RANSWallHeight",p.geometry.height,2.);
        session->LoadParameter("RANSWallRoughness",p.geometry.roughness,.1);session->LoadParameter("Kinvis",p.geometry.molecularNu,1.5e-5);
        session->LoadParameter("RANSWallSampleHeight",p.geometry.sampleHeight,p.geometry.height);session->LoadParameter("RANSWallTailStart",p.geometry.tailStart,p.geometry.height);
        session->LoadParameter("RANSWallNormalOrder",p.normalOrder,32);session->LoadParameter("RANSWallTangentialOrder",p.tangentialOrder,8);
        std::string model,policy,wall;session->LoadSolverInfo("TurbulenceModel",model,"MISSING");
        session->LoadSolverInfo("RANSKappaPolicy",policy,"Natural");session->LoadSolverInfo("RANSWallFormulation",wall,"MISSING");
        Require(policy=="Natural" || policy=="Prescribed","unknown kappa policy");
        Require(wall=="EPAWF" || wall=="AcceptedBudgetAWF","new-model wall formulation required");
        double kappa;session->LoadParameter("RANSKappa",kappa,.41);
        p.model=RANSResolvedModel::Resolve(ParseRANSModel(model),policy=="Natural"?RANSKappaPolicy::Natural:RANSKappaPolicy::Prescribed,
            kappa,wall=="EPAWF"?RANSWallFormulation::EPAWF:RANSWallFormulation::AcceptedBudgetAWF,p.geometry.molecularNu);
        Require(p.model.NewModel(),"new-model descriptor required; no CLS-default replay");
        for(const auto &name:{"BodyForceInner","BodyForceOuter","ScalarSourceInner","ScalarSourceOuter"})
            Require(!session->DefinesFunction(name),"neutral preservation verifier requires unforced source definitions");
        std::string referenceMode;session->LoadSolverInfo("RANSWallReference",referenceMode,"Neutral");
        Require(referenceMode=="Neutral","neutral verification requires the Neutral reference/boundary contract");
        double eta;session->LoadParameter("RANSWallEta",eta,0.);Require(eta==0,"unanchored verification requires eta zero");
        RANSWallEvolutionPatch patch(session,graph,world->GetComm(),p,LoadRANSWallNormalMap(session,2,p.geometry.height));
        Require(patch.Dimension()==2 && !patch.Mapped() && (patch.AcceptedAWF() || patch.UniformEPFlat()),
            "neutral preservation probe supports accepted or EP plateau references, serial flat2D only");
        RANSWallEvolutionData data(2);Vector state,reference;double time=request.get<double>("time",0.);
        const std::string operation=request.get<std::string>("operation","neutral");
        if(operation=="neutral")
        {
            double us,dx,dy;session->LoadParameter("RANSWallFrictionVelocity",us,.4);
            session->LoadParameter("RANSWallDirectionX",dx,1.);session->LoadParameter("RANSWallDirectionY",dy,0.);
            state=patch.NeutralReference(us,{dx,dy});data.fixedBoundary=patch.FixedReferenceBoundary(state);
            state=patch.Lift(state,time,data);reference=state;
        }
        else
        {
            Require(operation=="replay","unknown read-only verification operation");
            Require(request.get<std::string>("identity")==patch.PreservationIdentity(),"verification physics/mesh/space identity mismatch");
            Require(request.get<std::string>("descriptor")==p.model.Signature(),"verification resolved descriptor mismatch");
            state=Values(request,"state");reference=Values(request,"reference");
            Require(state.size()==std::size_t(patch.Size()) && reference.size()==state.size(),"verification coefficient size mismatch");
            for(const auto &item:request.get_child("fixed_boundary"))
            {auto q=item.second.begin();const int row=q->second.get_value<int>();++q;const double value=q->second.get_value<double>();Require(row>=0 && row<patch.Size(),"invalid fixed boundary row");Require(data.fixedBoundary.emplace(row,value).second,"duplicate fixed boundary row");}
            Require(data.fixedBoundary==patch.FixedReferenceBoundary(reference),"verification external boundary identity mismatch");
        }
        patch.SetReference(reference,data,time,0.);
        std::vector<RANSWallCoordinate> grid;
        if(auto points=request.get_child_optional("grid"))for(const auto &item:*points)
        {RANSWallCoordinate q{};int i=0;for(const auto &v:item.second){Require(i<3,"invalid verification coordinate");q[i++]=v.second.get_value<double>();}Require(i==3,"invalid verification coordinate");grid.push_back(q);}
        const auto stateCopy=state;
        const auto observation=patch.PreservationPacket(state,time,data,request.get<bool>("export_basis",false),grid);
        Require(state==stateCopy,"read-only verification modified coefficients");
        const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        Tree observed;std::istringstream input(observation);boost::property_tree::read_json(input,observed);
        const auto &audit=observed.get_child("audit_alpha_zero");
        const bool physical=audit.get_child("accepted_awf_scalar_domain").front().second.get_value<double>()==1 &&
            audit.get_child_optional("accepted_awf_hard_admissible") &&
            audit.get_child("accepted_awf_hard_admissible").front().second.get_value<double>()==1;
        std::string initialHistory;
        if(auto target=request.get_optional<std::string>("write_initial_history"))
        {
            Require(operation=="replay" && time==0. && physical &&
                std::memcmp(state.data(),reference.data(),state.size()*sizeof(double))==0,
                "initial history export requires an admissible identical state/reference at time zero");
            const auto scheme=session->GetTimeIntScheme();
            Require(scheme.method=="IMEX" && scheme.variant.empty() && scheme.order==1,
                "initial history export supports native one-sweep IMEX1 only");
            double dt,us,dx,dy;session->LoadParameter("TimeStep",dt);
            session->LoadParameter("RANSWallFrictionVelocity",us,.4);session->LoadParameter("RANSWallDirectionX",dx,1.);session->LoadParameter("RANSWallDirectionY",dy,0.);
            RANSWallTimeIntegration integrator(1,dt,
                [](const Vector&,double,double)->Vector{throw std::runtime_error("initial history export attempted a governing solve");},
                [&](const Vector &v,double t){auto full=state;std::copy(v.begin(),v.end(),full.begin());return patch.ExplicitBody(full,t,data);});
            integrator.Initialise(Vector(state.begin(),state.begin()+patch.DifferentialSize()),0.);
            RANSWallTimeIntegration::Checkpoint checkpoint;
            checkpoint.signature=RANSWallPhysicsSignature(patch,session,0.,true,us,{dx,dy},false);
            checkpoint.integration=integrator.Capture();checkpoint.coupledState=state;
            checkpoint.fixedReference=reference;checkpoint.fixedDefect=patch.Defect();checkpoint.eta=0.;
            initialHistory=RANSWallTimeIntegration::WriteCheckpoint(*target,checkpoint);
            Require(state==stateCopy,"initial history export modified coefficients");
        }
        const std::string output=J::Object({{"status",J::Quote(physical?"pass":"fail")},{"execution",J::Quote("completed")},
            {"physical_admissibility",J::Quote(physical?"pass":"fail")},
            {"accuracy",J::Quote("not_assessed_by_sampling_adapter")},{"performed_time_steps","0"},{"initial_history",initialHistory.empty()?"null":J::Quote(initialHistory)},
            {"reference",J::Numbers(reference)},{"fixed_boundary",J::Constraints(data.fixedBoundary)},
            {"elapsed_seconds",J::Number(elapsed)},{"observation",observation}});
        std::ofstream file(argv[3]);Require(bool(file),"cannot write verification output");file<<output<<'\n';
        return physical?0:1;
    }
    catch(const std::exception &error)
    {
        if(argc==4){std::ofstream file(argv[3]);file<<J::Object({{"status",J::Quote("fail")},{"failure",J::Quote(error.what())}})<<'\n';}
        std::cerr<<"preservation verification: "<<error.what()<<std::endl;return 1;
    }
}
