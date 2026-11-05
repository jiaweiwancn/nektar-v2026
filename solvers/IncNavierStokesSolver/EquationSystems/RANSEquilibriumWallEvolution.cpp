///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSEquilibriumWallEvolution.h>
#include <IncNavierStokesSolver/EquationSystems/IncNavierStokes.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <LibUtilities/Communication/CommSerial.h>
#include <LibUtilities/BasicUtils/FieldIO.h>
#include <SpatialDomains/MeshGraphIO.h>
#include <boost/property_tree/json_parser.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace Nektar
{
using namespace RANSWallAlgebra;
namespace
{
void Write(std::ostream &out,const Vector &v)
{out<<'[';for(int i=0;i<v.size();++i){if(i)out<<',';out<<v[i];}out<<']';}
void WriteRows(std::ostream &out,const std::vector<Vector> &v)
{out<<'[';for(int i=0;i<v.size();++i){if(i)out<<',';Write(out,v[i]);}out<<']';}
std::string Quoted(const std::string &value)
{
    std::ostringstream out;out<<'"';for(char c:value)
    {if(c=='"' || c=='\\')out<<'\\'<<c;else if(c=='\n')out<<"\\n";else out<<c;}out<<'"';return out.str();
}
}
std::string RANSEquilibriumWallEvolution::className=SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
    "RANSEquilibriumWallEvolution",RANSEquilibriumWallEvolution::create);
std::string RANSEquilibriumWallEvolution::solverTypeLookupId=LibUtilities::SessionReader::RegisterEnumValue(
    "SolverType","RANSEquilibriumWallEvolution",eRANSEquilibriumWallEvolution);
RANSWallEvolutionFunction RANSEquilibriumWallEvolution::Function(const std::string &name,const std::string &variable,bool required) const
{
    Require(!required || m_serialSession->DefinesFunction(name,variable),"missing evolution function "+name+"/"+variable);
    if(!m_serialSession->DefinesFunction(name,variable))return [](const RANSWallCoordinate&,double){return 0.;};
    Require(m_serialSession->GetFunctionType(name,variable)==LibUtilities::eFunctionTypeExpression,"evolution requires analytic session functions");
    const auto equation=m_serialSession->GetFunction(name,variable);const bool two=m_patch->Dimension()==2;
    return [equation,two,patch=m_patch,chart=m_chartFunctions](const RANSWallCoordinate &q,double t){
        const auto p=chart?q:patch->Physical(q);return equation->Evaluate(p[0],two?p[2]:p[1],two?0:p[2],t);};
}
std::vector<RANSWallEvolutionFunction> RANSEquilibriumWallEvolution::Functions(const std::string &name,bool required) const
{
    std::vector<RANSWallEvolutionFunction> result;const int dim=m_patch->Dimension();
    for(int b=0;b<m_patch->Gauge();++b)
    {
        const std::string variable=b<dim?m_session->GetVariable(b):(b==m_patch->K()?"k":(b==m_patch->Epsilon()?"eps_tilde":"p"));
        const auto inner=Function(name+"Inner",variable,required && b!=m_patch->Epsilon() && b!=m_patch->OuterPressure());
        const auto outer=Function(name+"Outer",variable,required && b!=m_patch->InnerPressure());
        if(b==m_patch->InnerPressure())result.push_back(inner);
        else if(b==m_patch->OuterPressure() || b==m_patch->Epsilon())result.push_back(outer);
        else result.push_back([inner,outer,H=m_patch->Parameters().geometry.height](const RANSWallCoordinate &p,double t){return p[2]<=H?inner(p,t):outer(p,t);});
    }
    return result;
}
void RANSEquilibriumWallEvolution::v_InitObject(bool declareFields)
{
    try
    {
        (void)declareFields;EquationSystem::v_InitObject(false);
        char program[]="nektar-wall-evolution";char *arguments[]{program,nullptr};
        // MeshGraphIO replaces GetFilenames() by each rank's partition. The
        // positional command-line inputs retain the complete original mesh.
        auto filenames=m_session->DefinesCmdLineArgument("input-file")?
            m_session->GetCmdLineArgument<std::vector<std::string>>("input-file"):m_session->GetFilenames();
        auto serial=LibUtilities::CommSerial::create(1,arguments);
        auto original=LibUtilities::SessionReader::CreateInstance(1,arguments,filenames,serial);original->InitSession();
        Require(!original->GetElement("NEKTAR/GEOMETRY")->Attribute("PARTITION"),
                "evolution requires the original unpartitioned XML input, not a partition directory");
        // Every collective rank uses the same global native CG numbering. Its
        // own source clone also keeps optimiser sidecars outside other ranks.
        const std::string nativeFile=m_sessionName+"-m4-native-rank"+std::to_string(m_session->GetComm()->GetRank())+".xml";
        TiXmlDocument document(original->GetDocument());auto root=document.FirstChildElement("NEKTAR");
        // Preserve effective world parameters, including command-line overrides,
        // in the full-mesh function interpreter. Keep their original case.
        if(auto parameters=root->FirstChildElement("CONDITIONS")->FirstChildElement("PARAMETERS"))
            for(auto item=parameters->FirstChildElement("P");item;item=item->NextSiblingElement("P"))
            {
                const std::string text=item->GetText()?item->GetText():"";const auto equals=text.find('=');
                if(equals==std::string::npos)continue;
                auto name=text.substr(0,equals);const auto first=name.find_first_not_of(" \t\r\n"),last=name.find_last_not_of(" \t\r\n");
                if(first==std::string::npos)continue;
                name=name.substr(first,last-first+1);
                if(m_session->DefinesParameter(name))
                {
                    std::ostringstream value;value<<std::setprecision(17)<<name<<" = "<<m_session->GetParameter(name);
                    item->Clear();item->LinkEndChild(new TiXmlText(value.str()));
                }
            }
        if(auto collection=root->FirstChildElement("COLLECTIONS"))root->RemoveChild(collection);
        auto collection=new TiXmlElement("COLLECTIONS");collection->SetAttribute("DEFAULT","NoCollection");root->LinkEndChild(collection);
        auto solverInfo=root->FirstChildElement("CONDITIONS")->FirstChildElement("SOLVERINFO");
        TiXmlElement *global=nullptr;
        for(auto item=solverInfo->FirstChildElement("I");item;item=item->NextSiblingElement("I"))
            if(item->Attribute("PROPERTY") && std::string(item->Attribute("PROPERTY"))=="GlobalSysSoln")global=item;
        if(!global){global=new TiXmlElement("I");global->SetAttribute("PROPERTY","GlobalSysSoln");solverInfo->LinkEndChild(global);}
        global->SetAttribute("VALUE","DirectFull");
        Require(document.SaveFile(nativeFile),"cannot write rank-local native session clone");
        filenames={nativeFile};m_serialSession=LibUtilities::SessionReader::CreateInstance(1,arguments,filenames,serial);m_serialSession->InitSession();
        m_serialSession->SetSolverInfo("GlobalSysSoln","DirectFull");
        auto graph=SpatialDomains::MeshGraphIO::Read(m_serialSession);const int dim=graph->GetMeshDimension();
        const std::vector<std::string> expected=dim==2?std::vector<std::string>{"u","v","k","eps_tilde","p"}:
            std::vector<std::string>{"u","v","w","k","eps_tilde","p"};
        Require(m_session->GetVariables()==expected,"evolution variables must be native velocity components,k,eps_tilde,p");
        RANSWallPatchParameters parameters;
        m_session->LoadParameter("RANSWallHeight",parameters.geometry.height,2.);
        m_session->LoadParameter("RANSWallRoughness",parameters.geometry.roughness,.1);
        m_session->LoadParameter("RANSWallSampleHeight",parameters.geometry.sampleHeight,1.);
        m_session->LoadParameter("RANSWallTailStart",parameters.geometry.tailStart,1.);
        m_session->LoadParameter("Kinvis",parameters.geometry.molecularNu,1.5e-5);
        m_session->LoadParameter("RANSWallChi",parameters.chi,1.);
        m_session->LoadParameter("RANSWallTangentialOrder",parameters.tangentialOrder,8);
        m_session->LoadParameter("RANSWallNormalOrder",parameters.normalOrder,16);
        std::string model;m_session->LoadSolverInfo("TurbulenceModel",model,"CLSHighRe");
        Require(model=="CLSHighRe" || model=="StandardKEpsilon","invalid evolution turbulence model");
        if(model=="StandardKEpsilon")parameters.model=SolverUtils::RANSNeutralReference::Standard();
        m_patch=std::make_shared<RANSWallEvolutionPatch>(m_serialSession,graph,m_session->GetComm(),parameters,
            LoadRANSWallNormalMap(m_session,dim,parameters.geometry.height));
        std::string coordinates;m_session->LoadSolverInfo("RANSWallFunctionCoordinates",coordinates,"Physical");
        Require(coordinates=="Physical" || coordinates=="Chart","unknown wall function coordinate convention");m_chartFunctions=coordinates=="Chart";
        std::string mode,condensation;m_session->LoadSolverInfo("RANSWallReference",mode,"Neutral");
        Require(mode=="Neutral" || mode=="Functions","RANSWallReference must be Neutral or Functions");m_neutral=mode=="Neutral";
        m_session->LoadSolverInfo("RANSWallCondense",condensation,"False");Require(condensation=="True" || condensation=="False","invalid condensation selector");m_condensed=condensation=="True";
        m_session->LoadParameter("RANSWallEta",m_eta,0.);m_session->LoadParameter("RANSWallFrictionVelocity",m_ustar,.4);
        m_session->LoadParameter("RANSWallDirectionX",m_direction[0],1.);m_session->LoadParameter("RANSWallDirectionY",m_direction[1],0.);
        m_session->LoadParameter("RANSWallNonlinearTolerance",m_tolerance,2e-12);m_session->LoadParameter("RANSWallMaxIterations",m_maximumIterations,20);
        int diagnostics=0;m_session->LoadParameter("RANSWallDiagnostics",diagnostics,0);
        Require(diagnostics==0 || diagnostics==1,"RANSWallDiagnostics must be 0 or 1");m_diagnostics=diagnostics;
        std::string rangeGuard;m_session->LoadSolverInfo("RANSWallPhysicalRangeGuard",rangeGuard,"Off");
        Require(rangeGuard=="Off" || rangeGuard=="R1LocalV1" || rangeGuard=="R1DiagnosticV1" || rangeGuard=="R2LocalV1","unknown physical-range guard");
        m_rangeGuardName=rangeGuard;m_r2RangeGuard=rangeGuard=="R2LocalV1";
        m_physicalRangeGuard=rangeGuard!="Off";m_diagnosticRangeContinuation=rangeGuard=="R1DiagnosticV1";
        if(m_physicalRangeGuard)
        {
            const auto robust=m_patch->RobustParameters();
            Require(m_diagnostics && m_eta==0 && m_session->GetComm()->GetSize()==1,
                "R1 applicability guard requires serial unanchored state diagnostics");
            Require(m_patch->VariantName()=="Regularized" && robust==Vector({.4,0,.2,.1,0,1,1}),
                "R1 applicability guard is declared only for beta=.4 GammaScale=1 and the unchanged linear capped source");
        }
        Require(m_eta==0 || (m_neutral && m_eta==1 && m_patch->Mapping().PlanarConstantHeight()),"eta=1 is restricted to isolated flat neutral verification");
        Require(!m_patch->TransitionEnabled() || m_eta==0,"every upstream transition/terrain run requires eta=0");
        Require(m_session->DefinesTimeIntScheme(),"evolution requires an explicit native IMEX time-integration declaration");
        const auto scheme=m_session->GetTimeIntScheme();
        Require(scheme.method=="IMEX" && scheme.variant.empty() && (scheme.order==1 || scheme.order==2),"evolution supports native IMEX orders 1/2, empty variant");m_order=scheme.order;
        m_session->LoadSolverInfo("RANSWallStageTest",m_stageTest,"Off");
        Require(m_stageTest=="Off" || m_stageTest=="MatchedR2" || m_stageTest=="FrozenScalarR2" ||
            m_stageTest=="ExactTargetR3" || m_stageTest=="JacobianBenchmarkR3" || m_stageTest=="DirectionExperimentR4" || m_stageTest=="DirectionExperimentR5","unknown native stage-test mode");
        const bool r3=m_stageTest=="ExactTargetR3" || m_stageTest=="JacobianBenchmarkR3" || m_stageTest=="DirectionExperimentR4" || m_stageTest=="DirectionExperimentR5";
        Require(m_stageTest!="DirectionExperimentR4" || m_patch->CoupledR4(),
            "DirectionExperimentR4 requires explicit CoupledR4 globalization");
        Require(!m_patch->CoupledR4() || m_stageTest=="ExactTargetR3" || m_stageTest=="DirectionExperimentR4",
            "CoupledR4 is restricted to the unforced alpha400 research stage");
        Require((m_stageTest!="DirectionExperimentR5" || m_patch->JointRestorationR5()) &&
            (!m_patch->JointRestorationR5() || m_stageTest=="ExactTargetR3" || m_stageTest=="DirectionExperimentR5"),
            "JointRestorationR5 requires the frozen exact target or its local comparison");
        if(m_stageTest!="Off" && !r3)
            Require(m_patch->NonlinearPolicyName()=="ScalarFeasibleR2" && m_eta==0 && m_neutral &&
                m_order==2 && m_steps==0 && m_timestep==.0025 && m_maximumIterations<=12 && m_tolerance==2e-12 &&
                m_session->GetComm()->GetSize()==1,"R2 stage test requires the frozen serial neutral IMEX2 startup input");
        if(r3)
        {
            m_session->LoadSolverInfo("RANSWallR3PredictorFile",m_r3PredictorFile,"None");
            Require(m_patch->NonlinearPolicyName()=="ScalarFeasibleR3" && m_patch->InnerStressName()=="SaturatedR2" &&
                m_eta==0 && m_neutral && m_order==2 && m_steps==0 && m_timestep==.0025 &&
                m_maximumIterations>0 && m_maximumIterations<=16 && m_tolerance==2e-12 &&
                m_session->GetComm()->GetSize()==1 && m_r3PredictorFile!="None",
                "R3 requires its explicit serial exact-target input and bound nonconverged predictor");
        }
        m_session->LoadSolverInfo("RANSWallRestart",m_restart,"None");m_session->LoadSolverInfo("RANSWallCheckpoint",m_checkpointBase,"None");
        Require(m_timestep>0,"evolution requires positive fixed TimeStep");m_data=RANSWallEvolutionData(dim);
        for(int a=0;a<dim;++a)
        {m_data.innerBody[a]=Function("BodyForceInner",expected[a]);m_data.outerBody[a]=Function("BodyForceOuter",expected[a]);}
        m_data.innerKSource=Function("ScalarSourceInner","k");m_data.outerKSource=Function("ScalarSourceOuter","k");m_data.outerEpsilonSource=Function("ScalarSourceOuter","eps_tilde");
        if(!m_neutral)
        {
            for(int b=0;b<dim;++b)
            {
                Function("BoundaryValuesInner",expected[b],true);
                Function("BoundaryValuesOuter",expected[b],true);
            }
            Function("BoundaryValuesInner","k",true);Function("BoundaryValuesOuter","k",true);
            Function("BoundaryValuesOuter","eps_tilde",true);
            const auto boundary=Functions("BoundaryValues",false);
            for(int b=0;b<m_patch->InnerPressure();++b)m_data.boundary[b]=boundary[b];
        }
        m_fields=Array<OneD,MultiRegions::ExpListSharedPtr>(expected.size());
        for(int a=0;a<dim;++a)m_fields[a]=MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(*m_patch->outerVelocity->Field());
        m_fields[dim]=m_patch->fields[m_patch->K()]->Field();m_fields[dim+1]=m_patch->fields[m_patch->Epsilon()]->Field();m_fields[dim+2]=m_patch->fields[m_patch->OuterPressure()]->Field();
    }
    catch(const std::exception &error){std::cerr<<"EP-AWF evolution initialisation: "<<error.what()<<std::endl;throw;}
}
std::string RANSEquilibriumWallEvolution::PhysicsSignature() const
{
    std::ostringstream s;s<<std::setprecision(17)<<m_patch->Signature()<<":eta="<<m_eta<<":neutral="<<m_neutral<<":us="<<m_ustar<<":direction="<<m_direction[0]<<','<<m_direction[1]<<":chart_functions="<<m_chartFunctions;
    if(!m_neutral && m_patch->CappedDestruction() && !m_patch->Mapped() && !m_patch->TransitionEnabled())
        s<<":capped_initial_trace=1";
    for(const auto &[name,value]:m_serialSession->GetParameters())
        if(name!="NUMSTEPS" && name!="IO_CHECKSTEPS" && name!="IO_INFOSTEPS" && name!="IO_TIMER_LEVEL" && name!="RANSWALLOWNERSHIPSHIFT")s<<'|'<<name<<'='<<value;
    const std::vector<std::string> names{"InitialConditionsInner","InitialConditionsOuter","InitialPerturbationInner","InitialPerturbationOuter","BoundaryValuesInner","BoundaryValuesOuter",
        "BodyForceInner","BodyForceOuter","ScalarSourceInner","ScalarSourceOuter","ExactSolutionInner","ExactSolutionOuter"};
    for(const auto &name:names)for(const auto &variable:m_session->GetVariables())if(m_serialSession->DefinesFunction(name,variable))
        s<<'|'<<name<<'/'<<variable<<'='<<m_serialSession->GetFunction(name,variable)->GetExpression();
    return s.str();
}
void RANSEquilibriumWallEvolution::v_DoInitialise(bool dump)
{
    try
    {
        (void)dump;m_state=m_neutral?m_patch->NeutralReference(m_ustar,m_direction):m_patch->Project(Functions("InitialConditions"),m_time);
        if(m_neutral)m_data.fixedBoundary=m_patch->FixedReferenceBoundary(m_state);
        else if(m_patch->Mapped() || m_patch->TransitionEnabled() || m_patch->CappedDestruction())
            m_state=m_patch->Lift(m_state,m_time,m_data);
        if(m_neutral && m_patch->TransitionEnabled())m_state=m_patch->Lift(m_state,m_time,m_data);
        m_patch->SetReference(m_state,m_data,m_time,m_eta);m_initial=m_state;
        const auto residualStart=std::chrono::steady_clock::now();
        m_rawNorms=m_patch->Norms(m_patch->Residual(m_state,m_state,0,m_time,m_data,false),m_patch->Constraints(m_state,m_time,m_data),m_ustar);
        const double rawResidualSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-residualStart).count();
        m_correctedNorms=m_patch->Norms(m_patch->Residual(m_state,m_state,0,m_time,m_data),m_patch->Constraints(m_state,m_time,m_data),m_ustar);
        m_implicitData=m_data;
        if(m_patch->TransitionEnabled())m_implicitData.sensorBody=m_data.outerBody;
        for(auto *body:{&m_implicitData.innerBody,&m_implicitData.outerBody})for(auto &f:*body)f=[](const RANSWallCoordinate&,double){return 0.;};
        m_signature=PhysicsSignature();
        m_integrator=std::make_unique<RANSWallTimeIntegration>(m_order,m_timestep,
            [this](const Vector &h,double t,double lambda){return Implicit(h,t,lambda);},
            [this](const Vector &s,double t){return Explicit(s,t);});
        if(m_restart!="None")
        {
            const auto checkpoint=RANSWallTimeIntegration::ReadCheckpoint(m_restart,m_signature,m_order,m_timestep,m_patch->Size(),m_patch->DifferentialSize());
            Require(checkpoint.eta==m_eta,"checkpoint eta differs");
            m_patch->RestoreReference(checkpoint.fixedReference,checkpoint.fixedDefect,checkpoint.eta,m_data,m_time);
            m_state=checkpoint.coupledState;m_initial=checkpoint.fixedReference;
            m_integrator->Restore(checkpoint.integration);m_time=m_integrator->Time();
        }
        else
        {
            // Perturb the state after freezing the neutral reference/defect.
            // This never changes the reference, wall model or boundary data.
            if(m_serialSession->DefinesFunction("InitialPerturbationInner") || m_serialSession->DefinesFunction("InitialPerturbationOuter"))
            {
                const auto perturbation=m_patch->Project(Functions("InitialPerturbation",false),m_time);
                for(int i=0;i<m_patch->DifferentialSize();++i)m_state[i]+=perturbation[i];
                m_state=m_patch->Lift(m_state,m_time,m_data);
            }
            m_integrator->Initialise(Vector(m_state.begin(),m_state.begin()+m_patch->DifferentialSize()),m_time);
        }
        m_historyReady=true;
        if(m_patch->NonlinearPolicyName()=="ScalarFeasibleR2" || m_patch->NonlinearPolicyName()=="ScalarFeasibleR3")
        {
            const auto domainStart=std::chrono::steady_clock::now();
            const auto domain=m_patch->ScalarDomain(m_state,m_time,m_data,true);
            m_initialEvaluationSeconds={rawResidualSeconds,
                std::chrono::duration<double>(std::chrono::steady_clock::now()-domainStart).count()};
            m_initialScalarDomain=RANSWallEvolutionPatch::ScalarDomainJSON(domain,true);
            Require(domain.admissible,"initial R2 scalar inventory is outside its evaluation domain: "+domain.failure);
        }
        if(m_diagnostics)RecordAudit(m_state,m_state,0,m_time);
        UpdateFields();
    }
    catch(const std::exception &error)
    {
        m_failure=error.what();
        if(m_diagnostics && m_session->GetComm()->GetRank()==0)
        {
            if(m_historyReady)WriteReport();
            else
            {
                // No native history exists yet. Preserve the original error
                // without calling Capture() on an uninitialised integrator.
                std::ofstream out(m_sessionName+".evolution.json");out<<std::setprecision(17)
                    <<"{\"schema\":1,\"initialization_complete\":false,\"failure\":"<<Quoted(m_failure)
                    <<",\"accepted_states\":[],\"rejected_state\":";Write(out,m_state);out<<"}\n";
            }
        }
        std::cerr<<"EP-AWF evolution state/history: "<<error.what()<<std::endl;throw;
    }
}
RANSEquilibriumWallEvolution::Vector RANSEquilibriumWallEvolution::Implicit(const Vector &history,double time,double lambda)
{
    auto fullHistory=m_state;std::copy(history.begin(),history.end(),fullHistory.begin());int iterations=0;
    const auto candidate=m_patch->Solve(m_state,fullHistory,1/lambda,time,m_implicitData,m_condensed,iterations,m_tolerance,m_maximumIterations);
    if(m_diagnostics)RecordAudit(candidate,fullHistory,1/lambda,time);
    m_state=candidate;
    const auto residual=m_patch->Residual(m_state,fullHistory,1/lambda,time,m_implicitData);
    const auto norms=m_patch->Norms(residual,m_patch->Constraints(m_state,time,m_implicitData),m_ustar);
    Vector record{time,lambda,1/lambda,static_cast<double>(iterations),Maximum(norms),m_patch->AnchorEnergyWork(m_state)};
    record.insert(record.end(),norms.begin(),norms.end());m_stageAudits.push_back(record);
    m_lastHistory=std::move(fullHistory);m_stageTime=time;m_stageAlpha=1/lambda;
    return Vector(m_state.begin(),m_state.begin()+m_patch->DifferentialSize());
}
void RANSEquilibriumWallEvolution::RecordAudit(const Vector &state,const Vector &history,double alpha,double time)
{
    RANSWallEvolutionPatch::Diagnostic audit;
    try{audit=m_r2RangeGuard?m_patch->CompleteRestrictedAudit(state,history,alpha,time,alpha?m_implicitData:m_data):
        m_patch->Audit(state,history,alpha,time,alpha?m_implicitData:m_data,m_ustar,m_physicalRangeGuard);}
    catch(const std::exception &error)
    {
        if(m_physicalRangeGuard)
        {
            m_rejectedState=state;
            m_rangeEvents.push_back({time,alpha,state,history,{},
                {std::string("range_audit_evaluation_failed:")+error.what()},true});
            throw std::runtime_error("declared R1 physical-range audit failed; complete candidate/history saved");
        }
        throw;
    }
    if(m_physicalRangeGuard)
    {
        std::vector<std::string> reasons;bool finitePhysicalExit=false;
        if(m_r2RangeGuard && audit.at("r2_scalar_domain")[0]!=1)
        {
            m_rejectedState=state;
            m_rangeEvents.push_back({time,alpha,state,history,audit,{"complete_scalar_domain_exit"},true});
            throw std::runtime_error("complete R2 scalar-domain exit before native history advance");
        }
        if(audit.at("admissible")[0]!=1)
        {reasons.push_back("sampled_scalar_or_sampler_inadmissibility");finitePhysicalExit=true;}
        const double K0=m_patch->Parameters().model.Q()*m_ustar*m_ustar;
        for(const std::string name:{"source_volume_inner","source_volume_outer","physical_sample_K","physical_sample_outer",
            "prescribed_source_inner","mass_volume","interface_inner","interface_outer","reconstruction_dense","reconstruction_native"})
        {
            const auto &row=audit.at(name);
            if(row[0] && row[3]<-1e-10*K0)
            {reasons.push_back("covariance:"+name);finitePhysicalExit=true;}
        }
        const auto &range=audit.at("outer_constitutive_range");
        // This is the explicit R1 restricted research box, not a claim that
        // the rejected wider M8 box has been repaired. No field is modified.
        if(range[0]>.001)reasons.push_back("outer_gradient_box_exit");
        if(range[1]<.99 || range[2]>1.01)reasons.push_back("outer_k_box_exit");
        if(range[3]<.98 || range[4]>1.02)reasons.push_back("outer_epsilon_box_exit");
        if(range[5]<.99 || range[6]>1.01)reasons.push_back("sampled_K_box_exit");
        const auto repair=audit.find("robust_source");
        if(repair!=audit.end() && repair->second[6]<.5)
        {reasons.push_back("lambda_applicability_exit");finitePhysicalExit=true;}
        if(!reasons.empty())
        {
            const bool stop=!m_diagnosticRangeContinuation || finitePhysicalExit;
            if(m_rangeEvents.empty() || stop)m_rangeEvents.push_back({time,alpha,state,history,audit,reasons,stop});
            if(stop)
            {
                m_rejectedState=state;
                throw std::runtime_error("declared R1 physical-range exit; candidate saved before advancing native history");
            }
        }
        if(m_r2RangeGuard && audit.at("r2_qualified")[0]!=1)
        {
            m_rejectedState=state;
            m_rangeEvents.push_back({time,alpha,state,history,audit,{"complete_R2_physical_range_exit"},true});
            throw std::runtime_error("complete R2 physical-range exit before native history advance");
        }
    }
    m_stateAudits.push_back(std::move(audit));m_auditTimes.push_back(time);
    if(m_stateAudits.back().at("admissible")[0]!=1)
    {m_rejectedState=state;throw std::runtime_error("accepted-stage candidate failed sampled admissibility; see state_audits");}
    m_acceptedStates.push_back(state);
}
RANSEquilibriumWallEvolution::Vector RANSEquilibriumWallEvolution::Explicit(const Vector &state,double time)
{
    auto full=m_state;std::copy(state.begin(),state.end(),full.begin());const auto derivative=m_patch->ExplicitBody(full,time,m_data);
    m_explicitAudits.push_back({time,Maximum(derivative)});return derivative;
}
void RANSEquilibriumWallEvolution::v_DoSolve()
{
    if(m_stageTest=="ExactTargetR3" || m_stageTest=="JacobianBenchmarkR3" || m_stageTest=="DirectionExperimentR4" || m_stageTest=="DirectionExperimentR5"){RunR3StageTest();return;}
    if(m_stageTest!="Off"){RunRestrictedStageTest();return;}
    try
    {
        for(int step=0;step<m_steps;++step)
        {
            const auto diff=m_integrator->Advance();std::copy(diff.begin(),diff.end(),m_state.begin());m_time=m_integrator->Time();
            const auto drift=m_patch->StateNorms(m_state,m_initial,m_ustar);
            Vector record{m_time,static_cast<double>(m_integrator->Steps()),Maximum(drift),m_patch->AnchorEnergyWork(m_state)};
            record.insert(record.end(),drift.begin(),drift.end());m_stepAudits.push_back(record);
            if(m_session->GetComm()->GetRank()==0)std::cout<<"EP-AWF IMEX"<<m_order<<" step "<<m_integrator->Steps()<<" time "<<m_time<<" stage scaled residual "<<m_stageAudits.back()[4]<<std::endl;
            if(m_checksteps>0 && m_integrator->Steps()%m_checksteps==0 && m_checkpointBase!="None")Checkpoint();
        }
        if(m_checkpointBase!="None")Checkpoint();
        if(m_patch->Mapped() || m_patch->TransitionEnabled() || m_diagnostics)
        {
            m_roughnessBudget=m_patch->RoughnessBudget(m_state,m_time,m_data);
            const auto history=m_lastHistory.empty()?m_state:m_lastHistory;
            m_energyBudget=m_patch->EnergyBudget(m_state,history,m_stageAlpha,m_time,m_implicitData);
            if(!m_lastHistory.empty())
            {
                m_finalStageLoads=m_patch->Residual(m_state,history,m_stageAlpha,m_stageTime,m_implicitData,false,false);
                for(int region=0;region<2;++region)m_regionalStageLoads.push_back(
                    m_patch->Residual(m_state,history,m_stageAlpha,m_stageTime,m_implicitData,false,false,region));
                m_finalStageConstraints=m_patch->Constraints(m_state,m_stageTime,m_implicitData);
            }
        }
        UpdateFields();
    }
    catch(const std::exception &error)
    {
        m_failure=error.what();
        if(m_diagnostics && m_session->GetComm()->GetRank()==0)WriteReport();
        std::cerr<<"EP-AWF native timestep: "<<error.what()<<std::endl;throw;
    }
}
RANSEquilibriumWallEvolution::Vector RANSEquilibriumWallEvolution::LoadR3Predictor(
    const Vector &history,double time)
{
    namespace J=RANSWallSolveJson;namespace D=RANSWallR3;
    boost::property_tree::ptree packet;boost::property_tree::read_json(m_r3PredictorFile,packet);
    const bool retained=m_patch->GroundSafeguard();
    const double originAlpha=retained?400.:3200.;
    const std::string role=(m_patch->ExperimentR9() || m_patch->ExperimentR10() || m_patch->ExperimentR11())?"retained_r8_terminal_iterate_only":
        (m_patch->ExperimentR8() || m_patch->PilotR8())?"retained_r7_terminal_iterate_only":
        (m_patch->InteriorProbeR7() || m_patch->AdaptiveR7())?"retained_r6_terminal_iterate_only":
        m_patch->JointRestorationR5()?"retained_r4_terminal_iterate_only":
        retained?"retained_r3_internal_iterate_only":"nonconverged_auxiliary_predictor_only";
    Require(packet.get<int>("schema")==1 && packet.get<std::string>("role")==role &&
        packet.get<std::string>("encoding")=="little-endian IEEE754 binary64" &&
        packet.get<double>("origin.alpha")==originAlpha && packet.get<double>("origin.time")==.0025 &&
        !packet.get<bool>("origin.converged") && packet.get<std::string>("origin.inner_stress_model")=="SaturatedR2" &&
        packet.get<double>("target.alpha")==400. && packet.get<double>("target.time")==time &&
        packet.get<double>("target.eta")==0. && packet.get<double>("target.tolerance")==m_tolerance &&
        !packet.get<bool>("target.history_advanced"),"R3 predictor changed its historical/target scope");
    auto pin=[&](const std::string &key)
    {
        const auto path=packet.get<std::string>(key+".path");
        Require(std::filesystem::path(path).is_absolute() && std::filesystem::is_regular_file(path) &&
            std::filesystem::file_size(path)==packet.get<std::uintmax_t>(key+".bytes") &&
            D::HashFile(path)==packet.get<std::string>(key+".sha256"),"R3 predictor file binding mismatch: "+key);
        return path;
    };
    const auto stateFile=pin("state"),historyFile=pin("history");
    const auto executable=pin("current_executable"),xml=pin("current_target_xml");
    Require(std::filesystem::exists("/proc/self/exe") && std::filesystem::equivalent(executable,"/proc/self/exe") &&
        D::HashFile("/proc/self/exe")==packet.get<std::string>("current_executable.sha256"),
        "R3 predictor belongs to a different current executable");
    const auto &inputs=m_session->GetFilenames();
    Require(inputs.size()==1 && std::filesystem::equivalent(xml,inputs[0]),"R3 predictor belongs to a different target XML");
    Require(packet.get<std::string>("current_native_signature_sha256")==D::Hash(m_signature) &&
        packet.get<std::string>("expected_h0_sha256")==D::Hash(history) &&
        packet.get<std::string>("history.sha256")==D::Hash(history),"R3 predictor changed native physics/mesh or original h0");
    auto read=[&](const std::string &file)
    {
        Require(std::filesystem::file_size(file)==history.size()*sizeof(double),"R3 predictor vector length differs");
        Vector values(history.size());std::ifstream in(file,std::ios::binary);
        in.read(reinterpret_cast<char*>(values.data()),values.size()*sizeof(double));
        Require(in.good() && std::all_of(values.begin(),values.end(),[](double value){return std::isfinite(value);}),
            "R3 predictor vector is incomplete or nonfinite");
        return values;
    };
    Require(D::Hash(read(historyFile))==D::Hash(history),"R3 predictor history payload differs from original h0");
    auto candidate=read(stateFile);const auto before=D::Hash(candidate);
    const auto domain=m_patch->ScalarDomain(candidate,time,m_implicitData,true);
    Require(domain.admissible,"R3 predictor is outside current scalar evaluation domain: "+domain.failure);
    const auto constraints=m_patch->Constraints(candidate,time,m_implicitData);
    const auto residual=m_patch->Residual(candidate,history,400.,time,m_implicitData);
    const auto norms=m_patch->Norms(residual,constraints,m_ustar);
    if(retained)
    {
        const auto residualFile=pin("origin.residual");
        Require(D::Hash(residual)==D::Hash(read(residualFile)),
            "Retained target state did not reproduce the exact original residual");
    }
    const auto audit=m_patch->CompleteRestrictedAudit(candidate,history,400.,time,m_implicitData,retained);
    Require(retained || (audit.count("r2_qualified") && audit.at("r2_qualified")[0]==1. && audit.count("r2_covariance") &&
        *std::min_element(audit.at("r2_covariance").begin(),audit.at("r2_covariance").end())>=0.),
        "R3 predictor exits current finite covariance/constitutive applicability sets");
    Require(before==D::Hash(candidate),"R3 predictor was silently lifted or changed during revalidation");
    std::vector<std::pair<std::string,std::string>> compactAudit;
    for(const auto &[name,value]:audit)if(name.rfind("r2_",0)==0)compactAudit.push_back({name,J::Numbers(value)});
    m_r3PredictorBinding=J::Object({{"manifest",D::Pin(m_r3PredictorFile)},
        {"state",D::Pin(stateFile)},{"history",D::Pin(historyFile)},
        {"current_executable",D::Pin(executable)},{"current_target_xml",D::Pin(xml)},
        {"current_native_signature_sha256",J::Quote(D::Hash(m_signature))},
        {"original_h0_sha256",J::Quote(D::Hash(history))},{"origin_converged","false"},
        {"origin_alpha",J::Number(originAlpha)},{"target_alpha","400"},{"target_time",J::Number(time)},
        {"predictor_is_physical_time_sample","false"},{"history_advanced","false"}});
    m_r3PredictorRevalidation=J::Object({{"scalar_domain",RANSWallEvolutionPatch::ScalarDomainJSON(domain)},
        {"current_target_scaled_block_norms",J::Numbers(norms)},{"current_target_constraints",J::Constraints(constraints)},
        {"current_target_residual_sha256",J::Quote(D::Hash(residual))},{"finite_set_audit",J::Object(compactAudit)},
        {"predictor_coefficients_unchanged","true"},{"historical_auxiliary_remains_nonconverged","true"}});
    return candidate;
}

void RANSEquilibriumWallEvolution::RunR3StageTest()
{
    namespace J=RANSWallSolveJson;
    Require(m_timestep==.0025 && m_time==0. && m_eta==0. && m_restart=="None" && m_checkpointBase=="None",
        "R3 exact target requires the unchanged unforced native startup and original history");
    for(const std::string name:{"BodyForceInner","BodyForceOuter","ScalarSourceInner","ScalarSourceOuter"})
        Require(!m_serialSession->DefinesFunction(name),"R3 exact target forbids prescribed body/scalar forcing");
    const Vector history=m_state;const auto native=m_integrator->Capture();const double targetTime=.0025;
    Require(native.steps==0,"R3 predictor import cannot follow a physical time step");
    m_stageTestInitialHistory=J::Object({{"state",J::Numbers(m_state)},{"time",J::Number(m_time)},
        {"steps",J::Number(native.steps)},{"history",J::Rows(native.history)},{"times",J::Numbers(native.times)}});
    auto candidate=LoadR3Predictor(history,targetTime);
    const bool benchmark=m_stageTest=="JacobianBenchmarkR3";
    const bool localR4=m_stageTest=="DirectionExperimentR4";
    const bool localR5=m_stageTest=="DirectionExperimentR5";
    Require((!(m_patch->DirectionReplayR6() || m_patch->CompareR7() || m_patch->InteriorProbeR7() || m_patch->ExperimentR8() || m_patch->ExperimentR9()) || localR5) &&
        (!(m_patch->AdaptiveR7() || m_patch->PilotR8()) || (!localR5 && !benchmark)),
        "R7 comparison and continuation require their designated stage modes");
    RANSWallEvolutionPatch::WallMatrix frozenMatrix;
    std::string frozenAssembly="[]";
    RANSWallEvolutionPatch::R8Ray ray;
    if(m_patch->PilotR8() ||
        (localR5 && !m_patch->CompareR7() && !m_patch->ExperimentR9() &&
        (m_patch->DirectionReplayR6() || m_patch->RefinedDirectionR6())))
    {
        namespace D=RANSWallR3;
        boost::property_tree::ptree packet;boost::property_tree::read_json(m_r3PredictorFile,packet);
        const auto &binding=packet.get_child("frozen_action");
        const auto path=binding.get<std::string>("matrix.path");
        const auto residual=m_patch->Residual(candidate,history,400.,targetTime,m_implicitData);
        Require(binding.get<std::string>("state_sha256")==D::Hash(candidate) &&
            binding.get<std::string>("history_sha256")==D::Hash(history) &&
            binding.get<std::string>("residual_sha256")==D::Hash(residual) &&
            std::filesystem::file_size(path)==candidate.size()*candidate.size()*sizeof(double) &&
            D::HashFile(path)==binding.get<std::string>("matrix.sha256"),
            "R6 replay matrix is not bound to this state, history and residual");
        std::ifstream binary(path,std::ios::binary);
        frozenMatrix.assign(candidate.size(),Vector(candidate.size()));
        for(auto &row:frozenMatrix)
            binary.read(reinterpret_cast<char*>(row.data()),row.size()*sizeof(double));
        Require(binary.good(),"R6 matrix payload is incomplete");
        frozenAssembly=binding.get<std::string>("assembly_json");
        if(m_patch->ExperimentR10() || m_patch->ExperimentR11())
        {
            auto vector=[&](const std::string &key)
            {
                const auto &p=packet.get_child(key);
                const auto file=p.get<std::string>("path");
                Require(std::filesystem::file_size(file)==candidate.size()*sizeof(double) &&
                    D::HashFile(file)==p.get<std::string>("sha256"),"R10 component binding mismatch");
                Vector value(candidate.size());std::ifstream in(file,std::ios::binary);
                in.read(reinterpret_cast<char*>(value.data()),value.size()*sizeof(double));
                Require(in.good(),"R10 component vector incomplete");return value;
            };
            if(m_patch->ExperimentR11())
            {
                for(const auto &entry:packet.get_child("r11_replay"))
                {
                    const auto &item=entry.second;
                    // Property-tree arrays have unnamed children, so bind this entry directly.
                    auto readPin=[&](const std::string &name)
                    {
                        const auto &pin=item.get_child(name);
                        const auto file=pin.get<std::string>("path");
                        Require(std::filesystem::file_size(file)==candidate.size()*sizeof(double) &&
                            D::HashFile(file)==pin.get<std::string>("sha256"),"R11 replay payload mismatch");
                        Vector value(candidate.size());std::ifstream in(file,std::ios::binary);
                        in.read(reinterpret_cast<char*>(value.data()),value.size()*sizeof(double));
                        Require(in.good(),"R11 replay vector incomplete");return value;
                    };
                    const auto direction=readPin("direction");
                    const double fraction=item.get<double>("fraction");
                    Vector actual;
                    if(item.get_child_optional("actual"))
                    {
                        auto expected=candidate;
                        for(int j=0;j<expected.size();++j)expected[j]+=fraction*direction[j];
                        Require(expected==readPin("state"),"R11 retrospective state leaves the bound ray");
                        actual=readPin("actual");
                    }
                    ray.replay.push_back({item.get<std::string>("name"),direction,actual,fraction,item.get<bool>("anchor")});
                }
            }
            ray.destruction=vector("r10_source.value");
            ray.ascentRegression=vector("r10_source.ascent");
            const auto &p=packet.get_child("r10_source.action");
            const auto file=p.get<std::string>("path");
            Require(std::filesystem::file_size(file)==candidate.size()*candidate.size()*sizeof(double) &&
                D::HashFile(file)==p.get<std::string>("sha256"),"R10 source action binding mismatch");
            ray.destructionAction.assign(candidate.size(),Vector(candidate.size()));
            std::ifstream in(file,std::ios::binary);
            for(auto &row:ray.destructionAction)
                in.read(reinterpret_cast<char*>(row.data()),row.size()*sizeof(double));
            Require(in.good(),"R10 source action incomplete");
        }
        if(m_patch->ExperimentR8() || m_patch->PilotR8())
        {
            auto readVector=[&](const boost::property_tree::ptree &pin)
            {
                const auto name=pin.get<std::string>("path");
                Require(std::filesystem::file_size(name)==candidate.size()*sizeof(double) &&
                    D::HashFile(name)==pin.get<std::string>("sha256"),"R8 retained vector binding differs");
                Vector value(candidate.size());
                std::ifstream stream(name,std::ios::binary);
                stream.read(reinterpret_cast<char*>(value.data()),value.size()*sizeof(double));
                Require(stream.good(),"R8 retained vector is incomplete");
                return value;
            };
            ray.step=readVector(packet.get_child("r8_ray.step"));
            ray.action=readVector(packet.get_child("r8_ray.action"));
            for(const auto &entry:packet.get_child("r8_ray.trials"))
            {
                const auto &trial=entry.second;
                const double f=trial.get<double>("fraction");
                auto predicted=candidate;
                for(int i=0;i<predicted.size();++i)predicted[i]+=f*ray.step[i];
                Require(predicted==readVector(trial.get_child("state")),"R8 retained trial leaves the J25 ray");
                Require(ray.residuals.emplace(f,readVector(trial.get_child("residual"))).second,
                    "Duplicate R8 retained trial fraction");
            }
        }
    }
    int updateLimit=(localR4 || localR5)?1:m_maximumIterations;
    if(m_patch->JointRestorationR5() && !localR5)
    {
        std::string requested;
        m_session->LoadSolverInfo("RANSWallR5Updates",requested,"");
        Require(!requested.empty() && requested.size()<=2 &&
            requested.find_first_not_of("0123456789")==std::string::npos,
            "R5 continuation requires an explicit numerical update budget");
        updateLimit=std::stoi(requested);
        Require(updateLimit>0 && updateLimit<=(m_patch->PilotR8()?12:32),"Stage update budget exceeds its bound");
    }
    bool solverReturned=false;int iterations=0;std::string failure;
    try
    {
        candidate=m_patch->SolveScalarFeasibleR3(candidate,history,400.,targetTime,m_implicitData,
            iterations,m_tolerance,updateLimit,benchmark,localR4,localR5,frozenMatrix,frozenAssembly,ray);solverReturned=true;
    }
    catch(const std::exception &error)
    {
        failure=std::string(J::ExceptionType(error))+": "+error.what();
        if(!m_patch->LastFeasibleState().empty())candidate=m_patch->LastFeasibleState();
    }
    const auto domain=m_patch->ScalarDomain(candidate,targetTime,m_implicitData,true);
    Vector residual,norms;std::map<int,double> constraints;RANSWallEvolutionPatch::Diagnostic audit;
    std::string residualFailure,auditFailure;
    if(domain.admissible)
    {
        try
        {
            constraints=m_patch->Constraints(candidate,targetTime,m_implicitData);
            residual=m_patch->Residual(candidate,history,400.,targetTime,m_implicitData);
            norms=m_patch->Norms(residual,constraints,m_ustar);
        }
        catch(const std::exception &error){residualFailure=std::string(J::ExceptionType(error))+": "+error.what();}
        try{audit=m_patch->CompleteRestrictedAudit(candidate,history,400.,targetTime,m_implicitData);}
        catch(const std::exception &error){auditFailure=std::string(J::ExceptionType(error))+": "+error.what();}
    }
    const bool coupled=!benchmark && !localR4 && !localR5 && solverReturned && residualFailure.empty() && !norms.empty() && Maximum(norms)<m_tolerance;
    const bool physical=coupled && auditFailure.empty() && audit.count("r2_qualified") && audit.at("r2_qualified")[0]==1.;
    std::vector<std::pair<std::string,std::string>> entries;
    for(const auto &[name,value]:audit)entries.push_back({name,J::Numbers(value)});
    m_stageTestLevels.push_back(J::Object({{"alpha","400"},{"time",J::Number(targetTime)},
        {"attempted",J::Boolean(!benchmark)},{"benchmark_only",J::Boolean(benchmark)},
        {"benchmark_completed",J::Boolean(benchmark && solverReturned)},{"converged",J::Boolean(coupled)},
        {"true_coupled_residual_converged",J::Boolean(coupled)},
        {"frozen_velocity_pressure_diagnostic","false"},{"iterations",J::Number(iterations)},
        {"failure",J::Quote(failure)},{"residual_verification_failure",J::Quote(residualFailure)},
        {"qualification_audit_failure",J::Quote(auditFailure)},
        {"state",J::Numbers(candidate)},{"history",J::Numbers(history)},
        {"constraints",J::Constraints(constraints)},{"residual",J::Numbers(residual)},
        {"scaled_block_norms",J::Numbers(norms)},
        {"scalar_domain",RANSWallEvolutionPatch::ScalarDomainJSON(domain,true)},
        {"physical_range_passed",J::Boolean(physical)},{"audit",J::Object(entries)},
        {"native_history_advanced","false"}}));
    m_targetConverged=coupled;m_targetPhysical=physical;
    if(!benchmark && !coupled && ((!localR4 && !localR5) || !solverReturned))m_failure="Exact target unresolved: "+failure;
    if(benchmark && !solverReturned)m_failure="R3 complete-J benchmark failed: "+failure;
    // The stage/predictor remains separate from original native state/history.
    WriteReport();
}

void RANSEquilibriumWallEvolution::RunRestrictedStageTest()
{
    namespace J=RANSWallSolveJson;
    const Vector history=m_state;auto guess=history;bool skipAuxiliary=false;
    const auto initialNative=m_integrator->Capture();
    m_stageTestInitialHistory=J::Object({{"state",J::Numbers(m_state)},
        {"time",J::Number(m_time)},{"steps",J::Number(initialNative.steps)},
        {"history",J::Rows(initialNative.history)},{"times",J::Numbers(initialNative.times)}});
    const double targetTime=m_time+m_timestep;
    const std::vector<double> levels=m_stageTest=="FrozenScalarR2"?std::vector<double>{400.}:
        std::vector<double>{3200.,1600.,800.,400.};
    Require(m_timestep==.0025 && m_eta==0 && m_restart=="None" && m_checkpointBase=="None",
        "bounded R2 stage test requires frozen unforced startup without restart/checkpoint");
    for(const std::string name:{"BodyForceInner","BodyForceOuter","ScalarSourceInner","ScalarSourceOuter"})
        Require(!m_serialSession->DefinesFunction(name),"R2 matched stage forbids manufactured/body/scalar forcing");
    for(double alpha:levels)
    {
        if(skipAuxiliary && alpha!=400.)
        {
            m_stageTestLevels.push_back(J::Object({{"alpha",J::Number(alpha)},
                {"attempted","false"},{"reason",J::Quote("remaining auxiliary levels skipped after a failed search")}}));
            continue;
        }
        auto candidate=guess;bool converged=false;int iterations=0;std::string failure;
        try
        {
            candidate=m_stageTest=="FrozenScalarR2"?
                m_patch->SolveFrozenScalars(guess,history,alpha,targetTime,m_implicitData,iterations,m_tolerance,m_maximumIterations):
                m_patch->Solve(guess,history,alpha,targetTime,m_implicitData,false,iterations,m_tolerance,m_maximumIterations);
            converged=true;
        }
        catch(const std::exception &error)
        {
            failure=std::string(J::ExceptionType(error))+": "+error.what();
            if(!m_patch->LastFeasibleState().empty())candidate=m_patch->LastFeasibleState();
            skipAuxiliary=true;
        }
        const auto domain=m_patch->ScalarDomain(candidate,targetTime,m_implicitData,true);
        Vector residual,norms;RANSWallEvolutionPatch::Diagnostic audit;
        std::map<int,double> constraints;std::string residualFailure,auditFailure;
        if(domain.admissible)
        {
            try
            {
                constraints=m_patch->Constraints(candidate,targetTime,m_implicitData);
                residual=m_patch->Residual(candidate,history,alpha,targetTime,m_implicitData);
                norms=m_patch->Norms(residual,constraints,m_ustar);
            }
            catch(const std::exception &error)
            {residualFailure=std::string(J::ExceptionType(error))+": "+error.what();}
            // Extra constitutive-reconstruction diagnostics can use more
            // normalization anchors than the residual. Preserve an obtained
            // coupled root even if this separate qualification audit fails.
            try{audit=m_patch->CompleteRestrictedAudit(candidate,history,alpha,targetTime,m_implicitData);}
            catch(const std::exception &error)
            {auditFailure=std::string(J::ExceptionType(error))+": "+error.what();}
        }
        const bool trueCoupled=converged && residualFailure.empty() && !norms.empty() && Maximum(norms)<m_tolerance;
        const bool physical=trueCoupled && auditFailure.empty() && audit.count("r2_qualified") && audit.at("r2_qualified")[0]==1;
        std::vector<std::pair<std::string,std::string>> entries;
        for(const auto &[name,value]:audit)entries.push_back({name,J::Numbers(value)});
        m_stageTestLevels.push_back(J::Object({{"alpha",J::Number(alpha)},{"time",J::Number(targetTime)},
            {"attempted","true"},{"converged",J::Boolean(converged)},
            {"true_coupled_residual_converged",J::Boolean(trueCoupled)},
            {"frozen_velocity_pressure_diagnostic",J::Boolean(m_stageTest=="FrozenScalarR2")},
            {"iterations",J::Number(iterations)},{"failure",J::Quote(failure)},
            {"residual_verification_failure",J::Quote(residualFailure)},{"qualification_audit_failure",J::Quote(auditFailure)},
            {"state",J::Numbers(candidate)},{"history",J::Numbers(history)},
            {"constraints",J::Constraints(constraints)},
            {"residual",J::Numbers(residual)},{"scaled_block_norms",J::Numbers(norms)},
            {"scalar_domain",RANSWallEvolutionPatch::ScalarDomainJSON(domain,true)},
            {"physical_range_passed",J::Boolean(physical)},{"audit",J::Object(entries)},
            {"native_history_advanced","false"}}));
        if(converged)guess=candidate;
        if(alpha==400.)
        {
            m_targetConverged=trueCoupled;m_targetPhysical=physical;
            if(!trueCoupled)m_failure="bounded target stage unresolved: "+failure;
        }
        std::cout<<"EP-AWF R2 stage test alpha "<<alpha<<" converged "<<trueCoupled
            <<" physical range "<<physical<<" history unchanged"<<std::endl;
    }
    // m_state, m_time, the IMEX history and accepted physical-stage arrays
    // deliberately remain at initialization. Auxiliary roots are not steps.
    WriteReport();
}
void RANSEquilibriumWallEvolution::Checkpoint()
{
    int bad=0;
    if(m_session->GetComm()->GetRank()==0)
    {
        try
        {
            RANSWallTimeIntegration::Checkpoint checkpoint;checkpoint.signature=m_signature;checkpoint.integration=m_integrator->Capture();
            checkpoint.coupledState=m_state;checkpoint.fixedReference=m_patch->Reference();checkpoint.fixedDefect=m_patch->Defect();checkpoint.eta=m_eta;
            m_lastCheckpoint=RANSWallTimeIntegration::WriteCheckpoint(m_checkpointBase,checkpoint);
        }
        catch(const std::exception &error){std::cerr<<error.what()<<std::endl;bad=1;}
    }
    m_session->GetComm()->Bcast(bad,0);Require(!bad,"cannot publish full-history checkpoint");m_session->GetComm()->Block();
}
void RANSEquilibriumWallEvolution::UpdateFields()
{
    const int dim=m_patch->Dimension();
    for(int a=0;a<dim;++a)
    {
        const auto values=m_patch->outerVelocity->Project([&](const RANSWallCoordinate &p){return m_patch->fields[a]->Point(p).Value(m_state,m_patch->offset[a]);});
        Array<OneD,NekDouble> global(values.size(),values.data());
        m_patch->outerVelocity->Field()->GetLocalToGlobalMap()->GlobalToLocal(global,m_fields[a]->UpdateCoeffs());
        m_fields[a]->BwdTrans(m_fields[a]->GetCoeffs(),m_fields[a]->UpdatePhys());m_fields[a]->SetPhysState(true);
    }
    for(int b=dim;b<dim+3;++b)
    {
        const int block=b==dim+2?m_patch->OuterPressure():b;
        m_patch->fields[block]->Store(m_state,m_patch->offset[block]);
    }
}
void RANSEquilibriumWallEvolution::v_GenerateSummary(SolverUtils::SummaryList &summary)
{
    SolverUtils::AddSummaryItem(summary,"Wall evolution","native IMEX, mapped multiregion CG");
    SolverUtils::AddSummaryItem(summary,"Physical dimension",m_patch->Dimension());SolverUtils::AddSummaryItem(summary,"Global unknowns",m_patch->Size());
    SolverUtils::AddSummaryItem(summary,"MPI backend","replicated global coefficients, owned-cell collective assembly");
    SolverUtils::AddSummaryItem(summary,"Correction weights","chi="+std::to_string(m_patch->Parameters().chi)+", eta="+std::to_string(m_eta));
}
void RANSEquilibriumWallEvolution::v_EvaluateExactSolution(unsigned int field,Array<OneD,NekDouble> &out,NekDouble time)
{
    const int nq=m_fields[field]->GetTotPoints();out=Array<OneD,NekDouble>(nq,0.);
    const auto exact=Functions("ExactSolution",false);Array<OneD,NekDouble> x(nq),y(nq),z(nq);m_fields[field]->GetCoords(x,y,z);
    const int block=field==m_patch->Dimension()+2?m_patch->OuterPressure():static_cast<int>(field);
    for(int i=0;i<nq;++i)
    {
        const RANSWallCoordinate point{x[i],m_patch->Dimension()==3?y[i]:0.,m_patch->Dimension()==3?z[i]:y[i]};
        out[i]=m_neutral?m_patch->fields[block]->Point(point).Value(m_patch->Reference(),m_patch->offset[block]):exact[block](point,time);
    }
}
void RANSEquilibriumWallEvolution::v_Output()
{
    if(m_session->GetComm()->GetRank()!=0)return;
    auto io=LibUtilities::FieldIO::CreateDefault(m_serialSession);
    auto save=[&](const std::string &name,const std::shared_ptr<RANSWallGlobalField> &native,const Vector &values)
    {
        native->Store(values);auto field=native->Field();auto definitions=field->GetFieldDefinitions();std::vector<std::vector<NekDouble>> data(definitions.size());
        for(int j=0;j<definitions.size();++j){definitions[j]->m_fields.push_back(name);field->AppendFieldData(definitions[j],data[j],field->UpdateCoeffs());}
        LibUtilities::FieldMetaDataMap meta;meta["Time"]=std::to_string(m_time);meta["PressureConvention"]="pi=p/rho+2k/3";
        if(m_patch->Mapped()){meta["CoordinateRepresentation"]="normal chart; physical geometry in VTU";meta["PhysicalMap"]=m_patch->Mapping().Signature();}
        io->Write(m_sessionName+"-"+name+(m_patch->Mapped()?".chart.fld":".fld"),definitions,data,meta,false);
        if(m_patch->Mapped())WriteMappedField(name,native,values);
    };
    for(int b=0;b<m_patch->Gauge();++b)
    {
        const auto native=b<m_patch->Dimension()?m_patch->outerVelocity:m_patch->fields[b];
        const std::string name=b<m_patch->Dimension()?m_session->GetVariable(b):(b==m_patch->K()?"k":(b==m_patch->Epsilon()?"eps_tilde":(b==m_patch->InnerPressure()?"pi_w":"p")));
        const auto values=native->Project([&](const RANSWallCoordinate &p){return m_patch->fields[b]->Point(p).Value(m_state,m_patch->offset[b]);});save(name,native,values);
    }
    for(int b=0;b<m_patch->Dimension();++b)
    {
        const auto values=m_patch->innerVelocity->Project([&](const RANSWallCoordinate &p){return m_patch->fields[b]->Point(p).Value(m_state,m_patch->offset[b]);});
        save(m_session->GetVariable(b)+"_w",m_patch->innerVelocity,values);
    }
    WriteReport();
}
void RANSEquilibriumWallEvolution::WriteMappedField(const std::string &name,const std::shared_ptr<RANSWallGlobalField> &field,const Vector &values) const
{
    field->Store(values);const auto native=field->Field();std::ofstream out(m_sessionName+"-"+name+".vtu");out<<std::setprecision(17);
    const int dimension=m_patch->Dimension(),points=native->GetTotPoints();int cells=0;
    for(int e=0;e<native->GetExpSize();++e){const auto exp=native->GetExp(e);int count=1;for(int a=0;a<dimension;++a)count*=exp->GetNumPoints(a)-1;cells+=count;}
    out<<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\"><UnstructuredGrid><Piece NumberOfPoints=\""<<points<<"\" NumberOfCells=\""<<cells<<"\"><Points><DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    Array<OneD,NekDouble> x(points),y(points),z(points);native->GetCoords(x,y,z);
    for(int j=0;j<points;++j){const auto p=m_patch->Physical({x[j],dimension==3?y[j]:0.,dimension==3?z[j]:y[j]});out<<p[0]<<' '<<p[1]<<' '<<p[2]<<'\n';}
    out<<"</DataArray></Points><PointData Scalars=\""<<name<<"\"><DataArray type=\"Float64\" Name=\""<<name<<"\" format=\"ascii\">\n";
    for(int j=0;j<points;++j)out<<native->GetPhys()[j]<<'\n';
    out<<"</DataArray></PointData><Cells><DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for(int e=0;e<native->GetExpSize();++e)
    {
        const auto exp=native->GetExp(e);const int nx=exp->GetNumPoints(0),ny=exp->GetNumPoints(1),nz=dimension==3?exp->GetNumPoints(2):1,start=native->GetPhys_Offset(e);
        for(int k=0;k<(dimension==3?nz-1:1);++k)for(int j=0;j<ny-1;++j)for(int i=0;i<nx-1;++i)
        {
            const int a=start+i+nx*j+nx*ny*k;out<<a<<' '<<a+1<<' '<<a+1+nx<<' '<<a+nx;
            if(dimension==3)out<<' '<<a+nx*ny<<' '<<a+nx*ny+1<<' '<<a+nx*ny+1+nx<<' '<<a+nx*ny+nx;
            out<<'\n';
        }
    }
    out<<"</DataArray><DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">";
    for(int j=1;j<=cells;++j)out<<j*(dimension==3?8:4)<<' ';
    out<<"</DataArray><DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">";
    for(int j=0;j<cells;++j)out<<(dimension==3?12:9)<<' ';
    out<<"</DataArray></Cells></Piece></UnstructuredGrid></VTKFile>\n";
    Require(static_cast<bool>(out),"cannot write physical mapped VTU field");
}
void RANSEquilibriumWallEvolution::WriteReport() const
{
    // No collective calls here: only rank zero writes already assembled state.
    std::ofstream out(m_sessionName+".evolution.json");out<<std::setprecision(17);
    const auto snapshot=m_integrator->Capture();
    out<<"{\"schema\":1,\"native_time_map\":true,\"dimension\":"<<m_patch->Dimension()<<",\"ranks\":"<<m_session->GetComm()->GetSize()
       <<",\"time\":"<<m_time<<",\"dt\":"<<m_timestep<<",\"steps\":"<<snapshot.steps<<",\"order\":"<<m_order
       <<",\"chi\":"<<m_patch->Parameters().chi<<",\"eta\":"<<m_eta<<",\"neutral_reference\":"<<(m_neutral?"true":"false")
       <<",\"condensed\":"<<(m_condensed?"true":"false")<<",\"checkpoint\":"<<Quoted(m_lastCheckpoint)<<",\"signature\":"<<Quoted(m_signature);
    out<<",\"inner_stress_model\":"<<Quoted(m_patch->InnerStressName())
       <<",\"nonlinear_policy\":"<<Quoted(m_patch->NonlinearPolicyName());
    if(!m_initialScalarDomain.empty())out<<",\"initial_scalar_domain\":"<<m_initialScalarDomain;
    if(!m_initialEvaluationSeconds.empty())
        out<<",\"initial_evaluation_seconds\":"<<RANSWallSolveJson::Object({
            {"unanchored_residual_constraints_norms",RANSWallSolveJson::Number(m_initialEvaluationSeconds[0])},
            {"scalar_domain",RANSWallSolveJson::Number(m_initialEvaluationSeconds[1])}});
    if(m_stageTest!="Off")
    {
        out<<",\"stage_test\":"<<Quoted(m_stageTest)<<",\"stage_test_target_time\":"<<m_time+m_timestep
           <<",\"stage_test_native_history_advanced\":false,\"stage_test_target_converged\":"<<(m_targetConverged?"true":"false")
           <<",\"stage_test_target_physical_range\":"<<(m_targetPhysical?"true":"false")
           <<",\"stage_test_initial_native_history\":"<<m_stageTestInitialHistory
           <<",\"stage_test_levels\":"<<RANSWallSolveJson::List(m_stageTestLevels);
        if(m_stageTest=="ExactTargetR3" || m_stageTest=="JacobianBenchmarkR3" || m_stageTest=="DirectionExperimentR4" || m_stageTest=="DirectionExperimentR5")
            out<<",\"r3_action_policy\":"<<Quoted(m_patch->ActionPolicyName())
               <<",\"r3_progress_file\":"<<Quoted(std::filesystem::absolute(m_patch->R3ProgressFile()).string())
               <<",\"r3_predictor_binding\":"<<m_r3PredictorBinding
               <<",\"r3_predictor_revalidation\":"<<m_r3PredictorRevalidation;
    }
    const auto &p=m_patch->Parameters();
    out<<",\"parameters\":";Write(out,{p.geometry.molecularNu,p.geometry.height,p.geometry.roughness,p.geometry.sampleHeight,
        p.geometry.tailStart,p.model.Q()*m_ustar*m_ustar,p.model.Q(),p.model.Kappa()});
    out<<",\"quadrature_orders\":";Write(out,{static_cast<double>(p.tangentialOrder),static_cast<double>(p.normalOrder)});
    out<<",\"ustar\":"<<m_ustar<<",\"direction\":";Write(out,{m_direction[0],m_direction[1]});
    out<<",\"cell_ownership\":[";int ownerCount=0;
    for(const auto &e:m_patch->fields[0]->Elements())
    {if(ownerCount++)out<<',';Write(out,{static_cast<double>(e.geometry),static_cast<double>(m_patch->Owner(e.geometry))});}out<<']';
    out<<",\"offsets\":[";for(int b=0;b<m_patch->offset.size();++b){if(b)out<<',';out<<m_patch->offset[b];}out<<']';
    out<<",\"state\":";Write(out,m_state);out<<",\"initial\":";Write(out,m_initial);
    out<<",\"fixed_reference\":";Write(out,m_patch->Reference());out<<",\"fixed_defect\":";Write(out,m_patch->Defect());
    out<<",\"raw_reference_loads\":";Write(out,m_patch->RawReferenceLoads());
    out<<",\"reference_constraints\":[";int constraintCount=0;
    for(auto [row,value]:m_patch->ReferenceConstraints())
    {if(constraintCount++)out<<',';Write(out,{static_cast<double>(row),value});}out<<']';
    out<<",\"raw_reference_norms\":";Write(out,m_rawNorms);out<<",\"corrected_reference_norms\":";Write(out,m_correctedNorms);
    out<<",\"defect_integrals\":";Write(out,m_patch->DefectIntegrals());out<<",\"anchor_energy_work\":"<<m_patch->AnchorEnergyWork(m_state);
    out<<",\"stage_audits\":";WriteRows(out,m_stageAudits);out<<",\"step_audits\":";WriteRows(out,m_stepAudits);out<<",\"explicit_audits\":";WriteRows(out,m_explicitAudits);
    if(m_diagnostics)
    {
        if(m_patch->SolveDiagnosticsEnabled())
            out<<",\"nonlinear_diagnostic_file\":"<<Quoted(m_patch->SolveDiagnosticFile())
               <<",\"nonlinear_diagnostic_written\":"<<(m_patch->SolveDiagnosticWritten()?"true":"false");
        if(m_physicalRangeGuard)
        {
            out<<",\"physical_range_guard\":"<<Quoted(m_rangeGuardName)
               <<",\"diagnostic_range_continuation\":"<<(m_diagnosticRangeContinuation?"true":"false")
               <<",\"physical_range_inventory\":{\"source_samplers_interfaces_mass_body\":true,\"native_GLL\":true,"
                 "\"inner_dense_points_per_axis\":\"2*p+7\",\"outer_dense_points_per_axis\":\"2*p+3\","
                 "\"covariance_tolerance_over_K0\":1e-10,\"lambda_minimum\":0.5,"
                 "\"restricted_gradient_radius\":0.001,\"restricted_k_ratio\":[0.99,1.01],\"restricted_epsilon_ratio\":[0.98,1.02]}"
               <<",\"physical_range_events\":[";
            for(int e=0;e<m_rangeEvents.size();++e)
            {
                if(e)out<<',';
                const auto &event=m_rangeEvents[e];
                out<<"{\"time\":"<<event.time<<",\"alpha\":"<<event.alpha
                   <<",\"stopped_before_history_advance\":"<<(event.stopped?"true":"false")<<",\"reasons\":[";
                for(int i=0;i<event.reasons.size();++i){if(i)out<<',';out<<Quoted(event.reasons[i]);}
                out<<"],\"state\":";Write(out,event.state);out<<",\"history\":";Write(out,event.history);
                out<<",\"audit\":{";int count=0;
                for(const auto &[name,value]:event.audit){if(count++)out<<',';out<<Quoted(name)<<':';Write(out,value);}
                out<<"}}";
            }
            out<<']';
        }
        out<<",\"diagnostics_schema\":1,\"failure\":"<<Quoted(m_failure)<<",\"rejected_state\":";Write(out,m_rejectedState);
        out<<",\"accepted_states\":";WriteRows(out,m_acceptedStates);out<<",\"audit_times\":";Write(out,m_auditTimes);
        out<<",\"state_audits\":[";
        for(int i=0;i<m_stateAudits.size();++i)
        {if(i)out<<',';out<<'{';int count=0;for(const auto &[name,value]:m_stateAudits[i])
            {if(count++)out<<',';out<<Quoted(name)<<':';Write(out,value);}out<<'}';}out<<']';
    }
    out<<",\"native_history\":";WriteRows(out,snapshot.history);out<<",\"native_history_times\":";Write(out,snapshot.times);
    out<<",\"last_stage_history\":";Write(out,m_lastHistory);out<<",\"last_stage_time\":"<<m_stageTime<<",\"last_stage_alpha\":"<<m_stageAlpha;
    const auto &mapping=m_patch->Mapping();
    out<<",\"mapped\":"<<(m_patch->Mapped()?"true":"false")<<",\"function_coordinates\":"<<Quoted(m_chartFunctions?"Chart":"Physical")
       <<",\"geometry_signature\":"<<Quoted(mapping.Signature())<<",\"map_parameters\":";
    Write(out,{static_cast<double>(mapping.kind),static_cast<double>(mapping.dimension),mapping.height,mapping.heightSlopeX,mapping.heightSlopeY,
        mapping.heightAmplitude,mapping.heightWaveX,mapping.heightWaveY,mapping.slopeX,mapping.slopeY,mapping.radius,mapping.amplitude,mapping.lengthX,mapping.lengthY});
    out<<",\"boundary_flux_correction\":"<<m_patch->BoundaryFluxCorrection()<<",\"roughness_budget\":";Write(out,m_roughnessBudget);
    out<<",\"energy_budget\":";Write(out,m_energyBudget);
    out<<",\"final_stage_loads\":";Write(out,m_finalStageLoads);
    out<<",\"regional_stage_loads\":";WriteRows(out,m_regionalStageLoads);
    out<<",\"final_stage_constraints\":[";constraintCount=0;
    for(auto [row,value]:m_finalStageConstraints){if(constraintCount++)out<<',';Write(out,{static_cast<double>(row),value});}out<<']';
    out<<",\"upstream_transition\":"<<(m_patch->TransitionEnabled()?"true":"false");
    out<<",\"model_variant\":"<<Quoted(m_patch->VariantName())<<",\"robust_parameters\":";
    Write(out,m_patch->RobustParameters());
    const auto &robustFlow=m_patch->RobustReferenceFlow(),&robustNormal=m_patch->RobustReferenceNormal();
    out<<",\"robust_reference_flow\":";Write(out,Vector(robustFlow.begin(),robustFlow.end()));
    out<<",\"robust_reference_normal\":";Write(out,Vector(robustNormal.begin(),robustNormal.end()));
    if(m_patch->TransitionEnabled())
    {
        out<<",\"epsilon_corner_policy\":\"dependent matching trace has priority over side-function corner values\"";
        const auto &t=m_patch->TransitionControls();const auto &c=t.state;const auto distances=m_patch->TransitionDistances();
        out<<",\"transition_controls\":";Write(out,{t.start,t.end,distances[0],distances[1],distances[2],
            c.on,c.off,c.marginOff,c.marginOn,c.reynoldsOff,c.reynoldsOn});
        out<<",\"transition_weights\":";Write(out,Vector(c.weights.begin(),c.weights.end()));
        out<<",\"transition_reference\":";Write(out,m_patch->TransitionReference());
        const auto trace=m_patch->DissipationTrace(m_state,m_time,m_data);out<<",\"transition_trace\":[";
        for(int part=0;part<2;++part){if(part)out<<',';out<<'[';int count=0;
            for(auto [id,value]:trace[part]){if(count++)out<<',';Write(out,{static_cast<double>(id),value});}out<<']';}out<<']';
        out<<",\"transition_samples\":[";int count=0;
        for(const auto &q:m_patch->Surface(0,2,p.geometry.height,1,8))
        {
            if(count++){out<<',';}
            const auto input=m_patch->TransitionInput(m_state,q.test.coordinate,m_time,m_data);
            const auto v=m_patch->Transition(m_state,q.test.coordinate,m_time,m_data);
            out<<"{\"chart\":";Write(out,Vector(q.test.coordinate.begin(),q.test.coordinate.end()));
            out<<",\"gates\":";Write(out,{v.chi,v.spatial,v.state,v.admissibility,v.marginGate,v.reynoldsGate,v.d2});
            out<<",\"input\":";Write(out,{input.coordinate,input.budgetMargin,input.reynolds});
            const auto &s=input.sensors;out<<",\"sensors\":[";Write(out,Vector(s.velocity.begin(),s.velocity.end()));out<<',';
            Write(out,Vector(s.velocityDifference.begin(),s.velocityDifference.end()));out<<',';Write(out,{s.normal});out<<',';
            Write(out,{s.k});out<<',';Write(out,{s.epsilon});out<<',';Write(out,Vector(s.pressure.begin(),s.pressure.end()));out<<',';
            Write(out,Vector(s.convection.begin(),s.convection.end()));out<<"]}";
        }out<<']';
    }
    out<<",\"fields\":[";
    for(int b=0;b<m_patch->Gauge();++b)
    {
        if(b)out<<',';
        const auto &field=m_patch->fields[b];const auto map=field->Field()->GetLocalToGlobalMap();
        out<<"{\"size\":"<<field->Size()<<",\"mesh\":[";
        for(int e=0;e<field->Elements().size();++e)
        {
            if(e)out<<',';
            const auto &el=field->Elements()[e];out<<"{\"geometry\":"<<el.geometry<<",\"degree\":"<<el.degree<<",\"lower\":";
            Write(out,Vector(el.lower.begin(),el.lower.end()));out<<",\"upper\":";Write(out,Vector(el.upper.begin(),el.upper.end()));out<<",\"ids\":[";
            const int nc=field->Field()->GetExp(e)->GetNcoeffs(),o=field->Field()->GetCoeff_Offset(e);
            for(int j=0;j<nc;++j){if(j)out<<',';out<<map->GetLocalToGlobalMap(o+j);}out<<"],\"sign\":[";
            for(int j=0;j<nc;++j){if(j)out<<',';out<<map->GetLocalToGlobalSign(o+j);}out<<"]}";
        }
        out<<"]}";
    }
    out<<"]}\n";Require(static_cast<bool>(out),"cannot write evolution report");
}
}
