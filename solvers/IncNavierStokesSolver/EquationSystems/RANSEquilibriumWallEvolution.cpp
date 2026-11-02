///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSEquilibriumWallEvolution.h>
#include <IncNavierStokesSolver/EquationSystems/IncNavierStokes.h>
#include <LibUtilities/Communication/CommSerial.h>
#include <LibUtilities/BasicUtils/FieldIO.h>
#include <SpatialDomains/MeshGraphIO.h>
#include <algorithm>
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
        Require(m_eta==0 || (m_neutral && m_eta==1 && m_patch->Mapping().PlanarConstantHeight()),"eta=1 is restricted to isolated flat neutral verification");
        Require(m_session->DefinesTimeIntScheme(),"evolution requires an explicit native IMEX time-integration declaration");
        const auto scheme=m_session->GetTimeIntScheme();
        Require(scheme.method=="IMEX" && scheme.variant.empty() && (scheme.order==1 || scheme.order==2),"evolution supports native IMEX orders 1/2, empty variant");m_order=scheme.order;
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
    for(const auto &[name,value]:m_serialSession->GetParameters())
        if(name!="NUMSTEPS" && name!="IO_CHECKSTEPS" && name!="IO_INFOSTEPS" && name!="IO_TIMER_LEVEL" && name!="RANSWALLOWNERSHIPSHIFT")s<<'|'<<name<<'='<<value;
    const std::vector<std::string> names{"InitialConditionsInner","InitialConditionsOuter","BoundaryValuesInner","BoundaryValuesOuter",
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
        else if(m_patch->Mapped())m_state=m_patch->Lift(m_state,m_time,m_data);
        m_patch->SetReference(m_state,m_data,m_time,m_eta);m_initial=m_state;
        m_rawNorms=m_patch->Norms(m_patch->Residual(m_state,m_state,0,m_time,m_data,false),m_patch->Constraints(m_state,m_time,m_data),m_ustar);
        m_correctedNorms=m_patch->Norms(m_patch->Residual(m_state,m_state,0,m_time,m_data),m_patch->Constraints(m_state,m_time,m_data),m_ustar);
        m_implicitData=m_data;
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
        else m_integrator->Initialise(Vector(m_state.begin(),m_state.begin()+m_patch->DifferentialSize()),m_time);
        UpdateFields();
    }
    catch(const std::exception &error){std::cerr<<"EP-AWF evolution state/history: "<<error.what()<<std::endl;throw;}
}
RANSEquilibriumWallEvolution::Vector RANSEquilibriumWallEvolution::Implicit(const Vector &history,double time,double lambda)
{
    auto fullHistory=m_state;std::copy(history.begin(),history.end(),fullHistory.begin());int iterations=0;
    m_state=m_patch->Solve(m_state,fullHistory,1/lambda,time,m_implicitData,m_condensed,iterations,m_tolerance,m_maximumIterations);
    const auto residual=m_patch->Residual(m_state,fullHistory,1/lambda,time,m_implicitData);
    const auto norms=m_patch->Norms(residual,m_patch->Constraints(m_state,time,m_implicitData),m_ustar);
    Vector record{time,lambda,1/lambda,static_cast<double>(iterations),Maximum(norms),m_patch->AnchorEnergyWork(m_state)};
    record.insert(record.end(),norms.begin(),norms.end());m_stageAudits.push_back(record);
    m_lastHistory=std::move(fullHistory);m_stageTime=time;m_stageAlpha=1/lambda;
    return Vector(m_state.begin(),m_state.begin()+m_patch->DifferentialSize());
}
RANSEquilibriumWallEvolution::Vector RANSEquilibriumWallEvolution::Explicit(const Vector &state,double time)
{
    auto full=m_state;std::copy(state.begin(),state.end(),full.begin());const auto derivative=m_patch->ExplicitBody(full,time,m_data);
    m_explicitAudits.push_back({time,Maximum(derivative)});return derivative;
}
void RANSEquilibriumWallEvolution::v_DoSolve()
{
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
        if(m_patch->Mapped())
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
    catch(const std::exception &error){std::cerr<<"EP-AWF native timestep: "<<error.what()<<std::endl;throw;}
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
    const auto &p=m_patch->Parameters();
    out<<",\"parameters\":";Write(out,{p.geometry.molecularNu,p.geometry.height,p.geometry.roughness,p.geometry.sampleHeight,
        p.geometry.tailStart,p.model.Q()*m_ustar*m_ustar,p.model.Q(),p.model.Kappa()});
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
