///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSEquilibriumWallSystem.h>
#include <IncNavierStokesSolver/EquationSystems/IncNavierStokes.h>
#include <fstream>
#include <iomanip>

namespace Nektar
{
using namespace RANSWallAlgebra;
using Layout=RANSWallFieldLayout;
std::string RANSEquilibriumWallSystem::className=SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
    "RANSEquilibriumWallSystem",RANSEquilibriumWallSystem::create);
std::string RANSEquilibriumWallSystem::solverTypeLookupId=LibUtilities::SessionReader::RegisterEnumValue(
    "SolverType","RANSEquilibriumWallSystem",eRANSEquilibriumWallSystem);
RANSWallTimeFunction RANSEquilibriumWallSystem::Function(const std::string &name,const std::string &variable,bool required) const
{
    Require(!required || m_session->DefinesFunction(name,variable),"missing "+name+"/"+variable);
    if(!m_session->DefinesFunction(name,variable))return [](double,double,double){return 0.;};
    Require(m_session->GetFunctionType(name,variable)==LibUtilities::eFunctionTypeExpression,
            "M3 stage inputs require analytic session expressions");
    const auto equation=m_session->GetFunction(name,variable);
    return [equation](double x,double r,double t){return equation->Evaluate(x,r,0,t);};
}
std::array<RANSWallTimeFunction,6> RANSEquilibriumWallSystem::Functions(const std::string &name,bool required) const
{
    std::array<RANSWallTimeFunction,6> f;
    const std::array<std::string,6> variables{"u","v","k","eps_tilde","p","p"};
    for(int b=0;b<6;++b)
    {
        const auto inner=Function(name+"Inner",variables[b],required && b!=Layout::Epsilon && b!=Layout::OuterPressure);
        const auto outer=Function(name+"Outer",variables[b],required && b!=Layout::InnerPressure);
        if(b==Layout::InnerPressure)f[b]=inner;
        else if(b==Layout::Epsilon || b==Layout::OuterPressure)f[b]=outer;
        else f[b]=[inner,outer,H=m_layout->Height()](double x,double r,double t){return r<=H?inner(x,r,t):outer(x,r,t);};
    }
    return f;
}
void RANSEquilibriumWallSystem::v_InitObject(bool declareFields)
{
    (void)declareFields;EquationSystem::v_InitObject(false);
    Require(m_spacedim==2 && m_session->GetVariables()==std::vector<std::string>({"u","v","k","eps_tilde","p"}),
            "M3 requires physical mesh (x,y)=(x,z), variables u,v,k,eps_tilde,p");
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
    Require(model=="CLSHighRe" || model=="StandardKEpsilon","unknown wall-patch turbulence model");
    if(model=="StandardKEpsilon")parameters.model=SolverUtils::RANSNeutralReference::Standard();
    m_session->LoadParameter("RANSWallStageAlpha",m_alpha,0.);
    m_session->LoadParameter("RANSWallNonlinearTolerance",m_tolerance,2e-12);
    m_session->LoadParameter("RANSWallMaxIterations",m_maximumIterations,15);
    std::string condensation;m_session->LoadSolverInfo("RANSWallCondense",condensation,"False");
    Require(condensation=="True" || condensation=="False","RANSWallCondense must be True or False");m_condensed=condensation=="True";
    std::string mode;m_session->LoadSolverInfo("RANSWallStageMode",mode,"Solve");
    Require(mode=="Solve" || mode=="Residual","RANSWallStageMode must be Solve or Residual");m_auditOnly=mode=="Residual";
    Require(m_steps==0,"M3 accepts a single supplied stage; native NumSteps/time integration belongs to M4");
    m_layout=std::make_shared<Layout>(m_session,m_graph,parameters.geometry.height);
    m_patch=std::make_shared<RANSWallPatch>(m_layout,parameters);
    m_fields[0]=m_layout->outerVelocity->Field();
    m_fields[1]=MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(*m_layout->outerVelocity->Field());
    m_fields[2]=m_layout->k->Field();m_fields[3]=m_layout->epsilon->Field();m_fields[4]=m_layout->outerPressure->Field();
    const std::array<std::string,4> boundaryVariables{"u","v","k","eps_tilde"};
    for(int b=0;b<4;++b)
    {
        const auto outer=Function("BoundaryValuesOuter",boundaryVariables[b],true);
        if(b==Layout::Epsilon)m_data.boundary[b]=outer;
        else
        {
            const auto inner=Function("BoundaryValuesInner",boundaryVariables[b],true);
            m_data.boundary[b]=[inner,outer,H=m_layout->Height()](double x,double r,double t){return r<=H?inner(x,r,t):outer(x,r,t);};
        }
    }
    for(int b=0;b<2;++b)
    {m_data.innerBody[b]=Function("BodyForceInner",b?"v":"u");m_data.outerBody[b]=Function("BodyForceOuter",b?"v":"u");}
    m_data.innerKSource=Function("ScalarSourceInner","k");
    m_data.outerKSource=Function("ScalarSourceOuter","k");
    m_data.outerEpsilonSource=Function("ScalarSourceOuter","eps_tilde");
}
void RANSEquilibriumWallSystem::v_DoInitialise(bool dump)
{
    (void)dump;m_state=m_patch->Project(Functions("InitialConditions",true),m_time);
    m_history=m_state;
    if(m_session->DefinesFunction("StageHistoryInner") || m_session->DefinesFunction("StageHistoryOuter"))
        m_history=m_patch->Project(Functions("StageHistory",true),m_time);
    UpdateFields();
}
void RANSEquilibriumWallSystem::v_DoSolve()
{
    if(!m_auditOnly)m_state=m_patch->Solve(m_state,m_history,m_alpha,m_time,m_data,m_condensed,m_iterations,m_tolerance,m_maximumIterations);
    m_residual=m_patch->Residual(m_state,m_history,m_alpha,m_time,m_data);
    m_norms=m_patch->ResidualNorms(m_residual,m_patch->Constraints(m_state,m_time,m_data));UpdateFields();
    std::cout<<(m_auditOnly?"EP-AWF flat residual audit, ":"EP-AWF flat stage converged, ")
             <<m_iterations<<" Newton iterations; max physical block dual norm "
             <<Maximum(m_norms)<<std::endl;
}
void RANSEquilibriumWallSystem::UpdateFields()
{
    m_patch->Store(m_state);
    auto values=m_layout->outerVelocity->Project([&](double x,double r){return m_layout->velocity->Point(x,r).Value(m_state,m_layout->offset[Layout::W]);});
    Array<OneD,NekDouble> global(values.size(),values.data());
    m_layout->outerVelocity->Field()->GetLocalToGlobalMap()->GlobalToLocal(global,m_fields[1]->UpdateCoeffs());
    m_fields[1]->BwdTrans(m_fields[1]->GetCoeffs(),m_fields[1]->UpdatePhys());m_fields[1]->SetPhysState(true);
}
void RANSEquilibriumWallSystem::v_GenerateSummary(SolverUtils::SummaryList &summary)
{
    SolverUtils::AddSummaryItem(summary,"Wall formulation","revised EP-AWF native flat mixed stage");
    SolverUtils::AddSummaryItem(summary,"Stage unknowns",m_layout->Size());
    SolverUtils::AddSummaryItem(summary,"Condensation",m_condensed?"inner free velocity; all pressure retained":"uncondensed");
    SolverUtils::AddSummaryItem(summary,"Pressure","pi=p/rho+2k/3, separate subdomain spaces, one volume mean");
}
void RANSEquilibriumWallSystem::v_EvaluateExactSolution(unsigned int field,Array<OneD,NekDouble> &output,NekDouble time)
{
    const int nq=m_fields[field]->GetTotPoints();output=Array<OneD,NekDouble>(nq,0.);
    if(!m_session->DefinesFunction("ExactSolutionInner") && !m_session->DefinesFunction("ExactSolutionOuter"))return;
    const auto exact=Functions("ExactSolution",true);
    Array<OneD,NekDouble> x(nq),r(nq),z(nq);m_fields[field]->GetCoords(x,r,z);
    const int block=field==4?Layout::OuterPressure:static_cast<int>(field);
    for(int q=0;q<nq;++q)output[q]=exact[block](x[q],r[q],time);
}
void RANSEquilibriumWallSystem::v_Output()
{
    // Each file keeps its own element set and coefficient count. Never expand
    // outer-only epsilon back onto the ground to fit the legacy output layout.
    for(int f=0;f<m_fields.size();++f)
    {
        std::vector<Array<OneD,NekDouble>> coefficients{m_fields[f]->UpdateCoeffs()};
        std::vector<std::string> variables{m_session->GetVariable(f)};
        auto field=m_fields[f];WriteFld(m_sessionName+"-"+variables[0]+".fld",field,coefficients,variables);
    }
    for(int component=0;component<3;++component)
    {
        auto native=component==2?m_layout->innerPressure:m_layout->innerVelocity;
        const int block=component==0?Layout::U:component==1?Layout::W:Layout::InnerPressure;
        auto values=native->Project([&](double x,double r){return m_layout->Field(static_cast<Layout::Block>(block)).Point(x,r,-1).Value(m_state,m_layout->offset[block]);});
        native->Store(values);auto field=std::static_pointer_cast<MultiRegions::ExpList>(native->Field());
        std::vector<Array<OneD,NekDouble>> coefficients{field->UpdateCoeffs()};
        std::vector<std::string> variables{component==0?"u_w":component==1?"v_w":"pi_w"};
        WriteFld(m_sessionName+"-"+variables[0]+".fld",field,coefficients,variables);
    }
    auto write=[](std::ostream &out,const Vector &v){out<<'[';for(int i=0;i<v.size();++i){if(i)out<<',';out<<v[i];}out<<']';};
    std::ofstream report(m_sessionName+".patch.json");report<<std::setprecision(17);
    report<<"{\"schema\":1,\"time\":"<<m_time<<",\"alpha\":"<<m_alpha<<",\"condensed\":"<<(m_condensed?"true":"false")
          <<",\"solved\":"<<(m_auditOnly?"false":"true")<<",\"iterations\":"<<m_iterations<<",\"native_time_map\":false,\"state\":";write(report,m_state);
    report<<",\"history\":";write(report,m_history);report<<",\"residual\":";write(report,m_residual);
    report<<",\"residual_norms\":";write(report,m_norms);
    report<<",\"history_as_state_residual\":";write(report,m_patch->Residual(m_history,m_history,m_alpha,m_time,m_data));
    report<<",\"forcing_samples\":[";
    int sampled=0;
    for(double fraction:{.13,.47,.91})
    {
        if(sampled++)report<<',';
        const double x=m_layout->Left()+fraction*(m_layout->Right()-m_layout->Left());
        write(report,{x,m_patch->SampleK(m_state,x),m_patch->Forcing(m_state,x,m_time,m_data)});
    }
    report<<']';
    report<<",\"source_probes\":[";
    for(int j=0;j<3;++j)
    {
        if(j)report<<',';
        const double x=(m_layout->Left()+m_layout->Right())/2,
            r=m_layout->Height()+(j+1)*(m_layout->Top()-m_layout->Height())/4;
        write(report,{x,r,m_data.outerKSource(x,r,m_time),m_data.outerEpsilonSource(x,r,m_time)});
    }
    report<<"]}\n";Require(static_cast<bool>(report),"cannot write wall-stage report");
}
}
