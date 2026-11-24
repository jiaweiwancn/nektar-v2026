// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeCoupled.h>
#include <IncNavierStokesSolver/EquationSystems/IncNavierStokes.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeWallV2.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <iostream>
#include <stdexcept>

#ifndef WJ_H27_SOURCE_ID
#error "The opt-in coupled executable requires a pinned scientific source identity"
#endif
namespace Nektar
{
namespace
{
double Clock(){return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();}
void Check(bool ok,const char*message){if(!ok)throw std::runtime_error(message);}
void VecJSON(std::ostream&out,const WJNative::Vec&v)
{out<<'[';for(int i=0;i<v.size();++i){if(i)out<<',';out<<RANSWallSolveJson::Number(v[i]);}out<<']';}
void AttemptJSON(std::ostream&out,const WJNative::CoupledStepInfo&info,const WJNative::CoupledOptions&options)
{
    out<<std::setprecision(17)<<"{\"nonlinear_iterations\":"<<info.nonlinearIterations<<",\"seconds\":"<<info.seconds<<",\"retries\":"<<info.retries
       <<",\"residuals\":["<<info.residual.momentum<<','<<info.residual.continuity<<','<<info.residual.gauge<<','<<info.residual.k<<','<<info.residual.epsilon<<','<<info.residual.trace
       <<"],\"coupled_krylov_iterations\":"<<info.coupledKrylovIterations<<",\"coupled_residual_evaluations\":"<<info.coupledResidualEvaluations
       <<",\"coupled_statistics_complete\":"<<(info.coupledStatisticsComplete?"true":"false")
       <<",\"coupled_preconditioner_applications\":"<<info.coupledPreconditionerApplications<<",\"rebuilds\":"<<info.rebuilds
       <<",\"disjoint_seconds\":{\"constitutive\":"<<info.timing.constitutive<<",\"wall\":"<<info.timing.wall
       <<",\"assembly\":"<<info.timing.assembly<<",\"flow\":"<<info.timing.flow<<",\"scalar\":"<<info.timing.scalar
       <<",\"rebuild\":"<<info.timing.rebuild<<",\"guards\":"<<info.timing.guards<<",\"forcing\":"<<info.timing.forcing
       <<",\"residual\":"<<info.timing.residual<<",\"nonlinear\":"<<info.timing.nonlinear<<",\"io\":"<<info.timing.io<<",\"other\":"<<info.timing.other<<'}'
       <<",\"acceleration_proposals\":"<<info.accelerationProposals<<",\"acceleration_accepted\":"<<info.accelerationAccepted<<",\"nonlinear_trace\":[";
    for(int i=0;i<int(info.nonlinearTrace.size());++i){const auto&t=info.nonlinearTrace[i];if(i)out<<',';
        out<<"{\"residuals\":["<<t.residual.momentum<<','<<t.residual.continuity<<','<<t.residual.gauge<<','<<t.residual.k<<','<<t.residual.epsilon<<','<<t.residual.trace
           <<"],\"accelerated\":"<<(t.accelerated?"true":"false")<<",\"seconds\":"<<t.seconds<<'}';}out<<']';
    WJNative::WriteStepPreconditionerFields(out,options,info);
    WJNative::WriteFlowSolverFields(out,options,info);out<<'}';
}
void SteadyJSON(std::ostream&out,const WJNative::SteadyKrylovResult&r,double mass)
{
    namespace J=RANSWallSolveJson;const auto&s=r.statistics;
    out<<std::setprecision(17)<<"{\"algorithm\":\""<<WJNative::SteadyKrylov::Identity<<"\",\"status\":"<<J::Quote(r.status)
       <<",\"reason\":"<<J::Quote(r.reason)<<",\"converged\":"<<(r.converged?"true":"false")
       <<",\"has_evaluated_state\":"<<(r.hasEvaluatedState?"true":"false")<<",\"true_norm\":"<<r.trueNorm
       <<",\"last_accepted_norm\":"<<r.lastAcceptedNorm<<",\"preconditioner_mass_only\":"<<mass
       <<",\"physical_time_advanced\":0,\"residual_evaluations\":"<<s.residualEvaluations<<",\"validity_evaluations\":"<<s.validityEvaluations
       <<",\"nonlinear_validation_enabled\":"<<(r.nonlinearValidationEnabled?"true":"false")<<",\"nonlinear_validation_evaluations\":"<<s.nonlinearValidationEvaluations<<",\"nonlinear_validation_rejected_trials\":"<<s.nonlinearValidationRejectedTrials
       <<",\"nonlinear_validation_seconds\":"<<s.nonlinearValidationSeconds
       <<",\"preconditioner_applications\":"<<s.preconditionerApplications<<",\"preconditioner_preparations\":"<<s.preconditionerPreparations
       <<",\"jacobian_actions\":"<<s.jacobianActions<<",\"fd_rejected_trials\":"<<s.fdRejectedTrials
       <<",\"line_search_rejected_trials\":"<<s.lineSearchRejectedTrials<<",\"linear_iterations\":"<<s.linearIterations
       <<",\"accepted_newton_steps\":"<<s.acceptedNewtonSteps<<",\"total_seconds\":"<<s.totalSeconds
       <<",\"callback_seconds\":{\"residual\":"<<s.residualSeconds<<",\"validity\":"<<s.validitySeconds<<",\"norm\":"<<s.normSeconds
       <<",\"precondition\":"<<s.preconditionerSeconds<<",\"prepare\":"<<s.preparationSeconds<<",\"arnoldi\":"<<s.arnoldiSeconds<<"},\"newton\":[";
    for(int i=0;i<int(r.newton.size());++i){const auto&n=r.newton[i];if(i)out<<',';out<<"{\"iteration\":"<<n.iteration<<",\"base_norm\":"<<n.baseNorm
        <<",\"forcing\":"<<n.forcing<<",\"linear_true_norm\":"<<n.linearTrueNorm<<",\"linear_converged\":"<<(n.linearConverged?"true":"false")
        <<",\"accepted\":"<<(n.accepted?"true":"false")<<",\"accepted_norm\":"<<n.acceptedNorm<<",\"accepted_step\":"<<n.acceptedStep<<'}';}
    out<<"],\"line_search\":[";for(int i=0;i<int(r.lineSearch.size());++i){const auto&l=r.lineSearch[i];if(i)out<<',';out<<"{\"newton\":"<<l.newton<<",\"trial\":"<<l.trial<<",\"step\":"<<l.step
        <<",\"nonlinear_validation_attempted\":"<<(l.nonlinearValidationAttempted?"true":"false")<<",\"nonlinear_validated\":"<<(l.nonlinearValidated?"true":"false")
        <<",\"admissible\":"<<(l.admissible?"true":"false")<<",\"evaluated\":"<<(l.evaluated?"true":"false")<<",\"true_norm\":"<<l.trueNorm
        <<",\"accepted\":"<<(l.accepted?"true":"false")<<",\"sufficient_decrease\":"<<(l.sufficientDecrease?"true":"false")<<",\"reason\":"<<J::Quote(l.reason)<<'}';}
    out<<"],\"best_from_accepted_iterate\":"<<(r.bestFromAcceptedIterate?"true":"false")<<",\"best_unknowns\":";VecJSON(out,r.x);out<<",\"best_raw_residual\":";VecJSON(out,r.residual);
    out<<",\"finite_differences\":[";for(int i=0;i<int(r.finiteDifferences.size());++i){const auto&f=r.finiteDifferences[i];if(i)out<<',';
        out<<"{\"newton\":"<<f.newton<<",\"cycle\":"<<f.cycle<<",\"iteration\":"<<f.iteration<<",\"reduction\":"<<f.reduction<<",\"purpose\":"<<J::Quote(f.purpose)<<",\"method\":"<<J::Quote(f.method)
           <<",\"relative_step\":"<<f.relativeStep<<",\"parameter_step\":"<<f.parameterStep<<",\"scaled_coordinate_step\":"<<f.scaledCoordinateStep<<",\"scaled_direction_infinity_norm\":"<<f.scaledDirectionInfinityNorm
           <<",\"nonlinear_validation_enabled\":"<<(r.nonlinearValidationEnabled?"true":"false")<<",\"admissibility_scope\":\"domain/root/covariance; uncommitted derivative probes\""
           <<",\"plus_admissible\":"<<(f.plusAdmissible?"true":"false")<<",\"minus_admissible\":"<<(f.minusAdmissible?"true":"false")<<",\"half_admissible\":"<<(f.halfAdmissible?"true":"false")
           <<",\"plus_reason\":"<<J::Quote(f.plusReason)<<",\"minus_reason\":"<<J::Quote(f.minusReason)<<",\"half_reason\":"<<J::Quote(f.halfReason)<<",\"seconds\":"<<f.seconds<<'}';}
    out<<"],\"linear_checks\":[";for(int i=0;i<int(r.linear.size());++i){const auto&l=r.linear[i];if(i)out<<',';
        out<<"{\"newton\":"<<l.newton<<",\"cycle\":"<<l.cycle<<",\"iteration\":"<<l.iteration<<",\"least_squares_rank\":"<<l.leastSquaresRank
           <<",\"projected_norm\":"<<l.projectedNorm<<",\"true_scaled_norm\":"<<l.trueScaledNorm<<",\"true_norm\":"<<l.trueNorm
           <<",\"true_checked\":"<<(l.trueChecked?"true":"false")<<",\"breakdown\":"<<(l.breakdown?"true":"false")<<",\"converged\":"<<(l.converged?"true":"false")<<'}';}out<<"]}";
}
}
std::string WJNativeCoupled::className=SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction("WJNativeCoupledV2",WJNativeCoupled::create);
std::string WJNativeCoupled::solverTypeLookupId=LibUtilities::SessionReader::RegisterEnumValue("SolverType","WJNativeCoupledV2",eWJNativeCoupledV2);
void WJNativeCoupled::v_InitObject(bool declareFields)
{
    try{
    const double start=Clock();(void)declareFields;EquationSystem::v_InitObject(false);
    Check(m_spacedim==2&&m_expdim==2&&m_HomogeneousType==eNotHomogeneous&&m_session->GetComm()->GetSize()==1,"WJNativeCoupledV2 requires serial 2D geometry");
    Check(m_session->GetVariables()==std::vector<std::string>({"u","v","k","epsilon","p"}),"WJNativeCoupledV2 variable order must be u,v,k,epsilon,p");
    Check(m_NumQuadPointsError==0,"restricted epsilon output currently requires NumQuadPointsError=0");
    std::string mode,table;m_session->LoadSolverInfo("WJNativeMode",mode,"Physical");Check(mode=="Physical"||mode=="Manufactured","unknown WJNativeMode");
    m_session->LoadSolverInfo("WJReferenceTable",table,"");Check(!table.empty(),"WJReferenceTable is required");
    std::string policy;m_session->LoadSolverInfo("WJReferencePolicy",policy,"Historical");
    auto reference=std::make_shared<WJNative::Reference>(table,WJNative::Reference::ParsePolicy(policy));
    double viscosity;m_session->LoadParameter("Kinvis",viscosity,reference->Parameters().nu);Check(std::isfinite(viscosity)&&std::abs(viscosity-reference->Parameters().nu)<=16*std::numeric_limits<double>::epsilon()*reference->Parameters().nu,"Kinvis differs from the pinned physical reference viscosity");
    WJNative::Geometry geometry;geometry.d0=reference->Parameters().d0;geometry.H=reference->Parameters().H;geometry.h=reference->Parameters().roughness;
    m_session->LoadParameter("WJPhysicalDisplacement",geometry.d0,geometry.d0);
    std::string geometryPolicy;
    m_session->LoadSolverInfo("WJGeometry",geometryPolicy,"Flat");
    Check(geometryPolicy=="Flat" || geometryPolicy=="NormalCosineHillV1","unknown WJGeometry");
    if (geometryPolicy=="NormalCosineHillV1")
    {
        auto &map=geometry.normalMap.emplace();
        map.kind=SolverUtils::RANSWallNormalMap::Kind::CosineHill;
        map.lengthX=2.5;
        map.displacement=geometry.d0;
        map.height=geometry.H;
        map.compositeTunnel=true;
        map.tunnelTop=mode=="Manufactured" ? geometry.d0+.5 : 20.;
        m_session->LoadParameter("WJHillHeight",map.amplitude,1.);
        Check(std::isfinite(map.amplitude) && map.amplitude>=0 && map.amplitude<=1.,
              "WJHillHeight must lie in [0,1] code units");
    }
    WJNative::CoupledOptions options;options.dt=m_timestep;options.manufactured=mode=="Manufactured";options.sourceIdentity=WJ_H27_SOURCE_ID;
    m_session->LoadParameter("WJHillInitialSlowdown",options.hillInitialSlowdown,0.);
    std::string backend;m_session->LoadSolverInfo("WJConstitutiveBackend",backend,"Historical");Check(backend=="Historical"||backend=="FastFull3DV1","unknown WJ constitutive backend");options.backend=backend=="Historical"?WJNative::ConstitutiveBackend::Historical:WJNative::ConstitutiveBackend::FastFull3DV1;
    m_session->LoadParameter("WJLinearTolerance",options.linearTolerance,1e-12);m_session->LoadParameter("WJCoupledTolerance",options.nonlinearTolerance,1e-10);m_session->LoadParameter("WJMaximumIterations",options.maximumIterations,60);
    m_session->LoadParameter("WJQuadratureX",options.quadratureX,8);m_session->LoadParameter("WJQuadratureR",options.quadratureR,12);
    m_session->LoadSolverInfo("WJStressSplit",options.stressSplit,"PositiveShearV1");
    m_session->LoadSolverInfo("WJSteadyVelocityPC",options.steadyVelocityPC,"FullGradientV1");
    m_session->LoadSolverInfo("WJImplicitSolver",options.advancement,"PicardV1");
    m_session->LoadSolverInfo("WJFlowPreconditioner", options.flowPreconditioner, "SparseLU");
    m_session->LoadParameter("WJPressureCorrectors", options.pressureCorrectors, 1);
    m_session->LoadSolverInfo("WJPicardStop", options.picardStop, "CoupledResidual");
    m_session->LoadSolverInfo("WJPicardAdvection", options.picardAdvection, "Explicit");
    m_session->LoadSolverInfo("WJNewtonAdvection", options.newtonAdvection, "Explicit");
    m_session->LoadSolverInfo("WJPicardStressLinearization", options.picardStressLinearization, "FrozenStrain");
    m_session->LoadParameter("WJPicardOuterCorrectors", options.picardOuterCorrectors, 0);
    m_session->LoadSolverInfo("WJStepPreconditioner",options.stepPreconditioner,"RebuildEveryNewtonV1");
    m_session->LoadParameter("WJStepPCRelativeDrift",options.stepPCRelativeDrift,.02);
    m_session->LoadParameter("WJStepPCMaximumUses",options.stepPCMaximumUses,8);
    if(options.advancement=="CoupledJFNKV1"){
        options.stepNewton.acceptanceNorm=options.nonlinearTolerance;
        m_session->LoadParameter("WJStepNewtonSeconds",options.stepNewton.maximumSeconds,0.);
        m_session->LoadParameter("WJStepNewtonIterations",options.stepNewton.maximumNewtonIterations,20);
        m_session->LoadParameter("WJStepKrylovDimension",options.stepNewton.restartDimension,40);
        m_session->LoadParameter("WJStepLinearCycles",options.stepNewton.maximumLinearCycles,8);
        m_session->LoadParameter("WJStepFDStep",options.stepNewton.fdRelativeStep,1e-6);
        double multiplier, ceiling;
        m_session->LoadParameter("WJStepForcingMultiplier", multiplier, 1.);
        m_session->LoadParameter("WJStepForcingCeiling", ceiling, .9);
        Check(std::isfinite(multiplier) && multiplier >= 1. && multiplier <= 4.,
              "WJStepForcingMultiplier must lie in [1,4]");
        Check(std::isfinite(ceiling) && ceiling >= .01 && ceiling <= .9,
              "WJStepForcingCeiling must lie in [.01,.9]");
        auto &newton = options.stepNewton;
        newton.forcingInitial = std::min(ceiling, multiplier * newton.forcingInitial);
        newton.forcingMinimum = std::min(ceiling, multiplier * newton.forcingMinimum);
        newton.forcingMaximum = std::min(ceiling, multiplier * newton.forcingMaximum);
        newton.forcingGamma *= multiplier;
    }
    std::string currentActions;m_session->LoadSolverInfo("WJCheckCurrentActions",currentActions,"False");Check(currentActions=="True"||currentActions=="False","invalid current action check switch");options.checkCurrentActions=currentActions=="True";
    m_session->LoadSolverInfo("WJOutletPolicy",options.outletPolicy,"FullTractionV1");
    m_session->LoadParameter("WJDiagnosticOutletX",options.diagnosticOutletX,15.);
    std::string outletActions;
    m_session->LoadSolverInfo("WJCheckOutletActions",outletActions,"False");
    Check(outletActions=="True"||outletActions=="False","invalid outlet action check switch");
    options.checkOutletActions=outletActions=="True";
    Check(!options.checkOutletActions||(!options.manufactured&&m_steps==0),
          "outlet action export requires the zero-time physical fixture");
    m_session->LoadParameter("WJAndersonDepth",options.andersonDepth,0);std::string trace;m_session->LoadSolverInfo("WJNonlinearTrace",trace,"False");Check(trace=="True"||trace=="False","WJNonlinearTrace must be True/False");options.recordNonlinearTrace=trace=="True";
    m_session->LoadParameter("WJMMSBodyAmplitude",options.mmsBodyAmplitude,0.);
    m_session->LoadSolverInfo("WJMMSVariant",options.mmsVariant,"Historical");
    m_session->LoadSolverInfo("WJInitialSourceIdentity",options.initialSourceIdentity,WJ_H27_SOURCE_ID);
    std::string solve;m_session->LoadSolverInfo("WJNativeSolve",solve,"Advance");Check(solve=="Advance"||solve=="SteadyReference","unknown WJ native solve mode");m_steadyMode=solve=="SteadyReference";
    Check(!m_steadyMode || options.flowPreconditioner == "SparseLU",
          "WJFlowPreconditioner=PressureCorrection is for time advancement");
    Check(!m_steadyMode || options.picardAdvection == "Explicit",
          "WJPicardAdvection=Implicit requires time advancement");
    Check(options.diagnosticOutletX==15.||m_steadyMode,"outlet location diagnostic requires steady preparation");
    if(m_steadyMode){Check(m_steps==0,"steady reference preparation requires NumSteps=0");
        m_session->LoadParameter("WJSteadySeconds",m_steadyOptions.maximumSeconds,0.);Check(m_steadyOptions.maximumSeconds>0,"steady reference preparation requires an explicit WJSteadySeconds cap");
        m_steadyOptions.acceptanceNorm=options.nonlinearTolerance;
        m_session->LoadParameter("WJSteadyNewtonIterations",m_steadyOptions.maximumNewtonIterations,50);
        m_session->LoadParameter("WJSteadyKrylovDimension",m_steadyOptions.restartDimension,40);
        m_session->LoadParameter("WJSteadyLinearCycles",m_steadyOptions.maximumLinearCycles,8);
        m_session->LoadParameter("WJSteadyFDStep",m_steadyOptions.fdRelativeStep,1e-6);
        m_session->LoadParameter("WJSteadyPreconditionerMass",m_steadyPreconditionerMass,0.);
        m_session->LoadParameter("WJSteadySeed",m_steadySeed,0.);Check(m_steadySeed>=0&&m_steadySeed<=.02,"invalid verification-only steady seed");Check(options.manufactured||m_steadySeed==0,"steady seed perturbation is verification-only");
    }
    m_session->LoadParameter("IO_InfoSteps",m_infoStride,20);m_session->LoadParameter("IO_CheckSteps",m_recoveryStride,200);m_session->LoadParameter("WJScientificSteps",m_scienceStride,200);m_session->LoadParameter("WJInjectRetryStep",m_injectRetryStep,-1);
    Check(m_infoStride>0&&m_recoveryStride>=0&&m_scienceStride>=0&&m_steps>=0,"invalid coupled output/time interval");
    m_session->LoadSolverInfo("WJOutputPrefix",m_prefix,m_sessionName+"-wj-v2");m_session->LoadSolverInfo("WJRestart",m_restart,"");std::string coefficients;m_session->LoadSolverInfo("WJOutputCoefficients",coefficients,"False");Check(coefficients=="True"||coefficients=="False","WJOutputCoefficients must be True/False");m_coefficients=coefficients=="True";
    std::string terminal;m_session->LoadSolverInfo("WJTerminalOutput",terminal,"False");Check(terminal=="True"||terminal=="False","WJTerminalOutput must be True/False");m_terminalOutput=terminal=="True";
    m_session->LoadParameter("WJTerminalBegin",m_terminalBegin,0.);m_session->LoadParameter("WJTerminalEnd",m_terminalEnd,0.);
    if(m_terminalOutput)Check(!m_steadyMode&&m_coefficients&&std::isfinite(m_terminalBegin)&&std::isfinite(m_terminalEnd)&&m_terminalBegin>=0&&m_terminalEnd>m_terminalBegin,"terminal output requires advancement, complete coefficients and a finite increasing nonnegative interval");
    m_session->LoadSolverInfo("WJInitialFields",m_initialFields,"");Check(m_initialFields.empty()||m_restart.empty(),"initial fields and cold restart are separate mutually exclusive inputs");
    m_session->LoadSolverInfo("WJPerturbationCertificate",m_perturbationCertificate,"");
    m_session->LoadSolverInfo("WJPerturbationCertificateSHA256",m_perturbationSha256,"");
    Check(m_perturbationCertificate.empty()==m_perturbationSha256.empty()&&
          (m_perturbationCertificate.empty()||(!options.manufactured&&!m_steadyMode&&
                                               !m_initialFields.empty()&&m_restart.empty())),
          "physical perturbation requires a certified initial-field import and fresh histories");
    std::vector<MultiRegions::ContFieldSharedPtr>fields(5);
    for(int i:{0,1,2,4}){fields[i]=MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(m_session,m_graph,m_session->GetVariable(i),true,false);m_fields[i]=fields[i];}
    std::vector<unsigned int>outer;
    for(int i=0;i<fields[2]->GetExpSize();++i){auto exp=fields[2]->GetExp(i);double low=1e100,high=-1e100;
        for(int j=0;j<exp->GetGeom()->GetNumVerts();++j){double x,z,d;exp->GetGeom()->GetVertex(j)->GetCoords(x,z,d);low=std::min(low,z-geometry.d0);high=std::max(high,z-geometry.d0);}
        Check(!(low<geometry.H-1e-12&&high>geometry.H+1e-12),"H must coincide with a complete native mesh edge layer");if(low>=geometry.H-1e-12)outer.push_back(i);
    }
    Check(!outer.empty(),"empty outer epsilon domain");auto subset=MemoryManager<MultiRegions::ExpList>::AllocateSharedPtr(*fields[2],outer,true);fields[3]=MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(m_session,*subset);m_fields[3]=fields[3];
    m_coupling=std::make_unique<WJNative::Coupling>(fields,geometry,reference,options);
    m_constructionSeconds=Clock()-start;
    }catch(const std::exception&e){std::cerr<<"WJ-H27 construction: "<<e.what()<<std::endl;throw std::runtime_error(e.what());}
}
void WJNativeCoupled::v_DoInitialise(bool dump)
{
    (void)dump;const double start=Clock();
    try{
        Check(!std::filesystem::exists(m_prefix+".steps.jsonl")&&!std::filesystem::exists(m_prefix+".result.json")&&!std::filesystem::exists(m_prefix+".science-index.json"),"coupled output prefix already exists");
        if (!m_restart.empty())
        {
            m_history = m_coupling->ReadHistory(m_restart, &m_recoveryTiming);
        }
        else if(!m_initialFields.empty())m_history=m_coupling->ReadInitialFields(m_initialFields,m_steadyMode);
        else m_history=m_coupling->Initialise();
        if(!m_perturbationCertificate.empty())m_history=m_coupling->PerturbInitialHistory(
            m_history,m_perturbationCertificate,m_perturbationSha256);
        if(m_steadySeed){Check(m_restart.empty(),"steady seed cannot modify a restart");auto seed=m_history.current;const auto&p=m_coupling->ReferenceModel().Parameters();
            for(int id:m_coupling->KKernel().FreeIds())seed.k[id]+=m_steadySeed*p.KIn*std::sin(.7*(id+1));
            for(int id:m_coupling->EpsilonKernel().FreeIds())seed.epsilon[id]+=m_steadySeed*m_coupling->ReferenceModel().Amplitude(p.KIn)/p.H*std::cos(.9*(id+1));
            m_history=m_coupling->FreshHistory(seed);
        }
        m_time=m_history.time;m_initialStep=m_history.step;m_initialTime=m_history.time;
        if (m_steps==0 && !m_steadyMode && m_history.step==0 &&
            m_coupling->FlowKernel().Core().Mapping().normalMap)
            m_last=m_coupling->AuditSteadyState(m_history.current);
        m_coupling->Store(m_history.current);m_stepsLog.open(m_prefix+".steps.jsonl");m_stepsLog<<std::setprecision(17);Check(bool(m_stepsLog),"cannot create coupled step log");
        m_initialisationSeconds=Clock()-start;
    }catch(const std::exception&e){std::cerr<<"WJ-H27 initialisation: "<<e.what()<<std::endl;throw std::runtime_error(e.what());}
}
unsigned WJNativeCoupled::TerminalReasons(double left,double right,bool before)const
{
    if(!m_terminalOutput)return 0;
    unsigned reasons=0;
    const std::array<double,2>boundary{{m_terminalBegin,m_terminalEnd}};
    for(int i=0;i<2;++i)
    {
        const double t=boundary[i];
        const unsigned lo=i?EndLeft:BeginLeft,hi=i?EndRight:BeginRight,exact=i?EndExact:BeginExact;
        if(before){if(t==left)reasons|=exact;else if(left<t&&t<right)reasons|=lo;}
        else{if(t==right)reasons|=exact;else if(left<t&&t<right)reasons|=hi;}
    }
    return reasons;
}
void WJNativeCoupled::WriteScientificFrame(unsigned reasons)
{
    if(!reasons)return;
    if(!m_scientificFrames.empty())
    {
        auto&last=m_scientificFrames.back();
        Check(last.step<=m_history.step,"scientific output clock moved backwards");
        if(last.step==m_history.step){Check(last.time==m_history.time,"duplicate scientific step has a different time");last.reasons|=reasons;return;}
    }
    const double began=Clock();ScientificFrame frame;frame.step=m_history.step;frame.time=m_history.time;frame.reasons=reasons;
    frame.path=m_prefix+".science-"+std::to_string(m_history.step)+".json";
    frame.origin=m_currentProcessAccepted?"current_process_accepted_step":m_restart.empty()?"fresh_initial_observation":m_history.step?"imported_accepted_history":"imported_initial_history";
    Check(!std::filesystem::exists(frame.path),"scientific output is immutable");
    std::ofstream out(frame.path);Check(bool(out),"cannot create scientific output");m_coupling->Report(out,m_history,m_last,m_coefficients);out<<'\n'<<std::flush;
    Check(bool(out),"scientific output write failed");out.close();Check(!out.fail(),"scientific output close failed");
    std::ifstream in(frame.path,std::ios::binary);Check(bool(in),"cannot reopen scientific output for identity");
    RANSWallR3::Sha256 digest;std::array<char,65536>buffer;
    while(in){in.read(buffer.data(),buffer.size());const auto count=in.gcount();if(count){digest.Update(buffer.data(),count);frame.bytes+=std::uint64_t(count);}}
    Check(in.eof()&&!in.bad(),"scientific output identity read failed");frame.sha256=digest.Finish();
    frame.seconds=Clock()-began;m_scientificFrameSeconds+=frame.seconds;m_scientificFrames.push_back(std::move(frame));
}
void WJNativeCoupled::WriteScientificIndex()
{
    namespace J=RANSWallSolveJson;const double began=Clock();const std::string path=m_prefix+".science-index.json";
    Check(!std::filesystem::exists(path),"scientific index is immutable");std::ofstream out(path);Check(bool(out),"cannot create scientific index");
    out<<std::setprecision(17)<<"{\"schema\":\"WJH27NativeScientificIndexV1\",\"source_identity\":"<<J::Quote(m_coupling->Options().sourceIdentity)
       <<",\"reference_identity\":"<<J::Quote(m_coupling->ReferenceModel().Identity())<<",\"signature\":"<<J::Quote(m_coupling->Signature())
       <<",\"run_complete\":"<<(m_solveCompleted?"true":"false")<<",\"dt\":"<<m_history.dt<<",\"initial_step\":"<<m_initialStep<<",\"initial_time\":"<<m_initialTime<<",\"final_step\":"<<m_history.step<<",\"final_time\":"<<m_history.time
       <<",\"coefficients_included\":"<<(m_coefficients?"true":"false")<<",\"terminal_output\":"<<(m_terminalOutput?"true":"false")<<",\"terminal_begin\":"<<m_terminalBegin<<",\"terminal_end\":"<<m_terminalEnd
       <<",\"reason_bits\":{\"periodic\":1,\"begin_left\":2,\"begin_right\":4,\"begin_exact\":8,\"end_left\":16,\"end_right\":32,\"end_exact\":64,\"final\":128},\"frames\":[";
    for(int i=0;i<int(m_scientificFrames.size());++i){const auto&f=m_scientificFrames[i];if(i)out<<',';
        out<<"{\"path\":"<<J::Quote(f.path)<<",\"sha256\":"<<J::Quote(f.sha256)<<",\"bytes\":"<<f.bytes<<",\"step\":"<<f.step<<",\"time\":"<<f.time
           <<",\"reason_mask\":"<<f.reasons<<",\"origin\":"<<J::Quote(f.origin)<<",\"write_seconds\":"<<f.seconds<<'}';}
    out<<"]}\n"<<std::flush;Check(bool(out),"scientific index write failed");out.close();Check(!out.fail(),"scientific index close failed");m_scientificIndexSeconds=Clock()-began;
}
void WJNativeCoupled::v_DoSolve()
{
    const double solveStart=Clock();
    try{
        if(m_steadyMode){
            m_steadyResult=m_coupling->SolveSteady(m_history,m_steadyOptions,m_last,m_steadyPreconditionerMass);
            const std::string path=m_prefix+".steady.json";Check(!std::filesystem::exists(path),"steady evidence is immutable");std::ofstream out(path);
            Check(bool(out),"cannot create steady evidence");SteadyJSON(out,m_steadyResult,m_steadyPreconditionerMass);out<<'\n'<<std::flush;Check(bool(out),"steady evidence write failed");out.close();Check(!out.fail(),"steady evidence close failed");
            if(!m_steadyResult.converged&&m_steadyResult.hasEvaluatedState){
                auto candidate=m_coupling->SteadyCandidateHistory(m_steadyResult.x);auto diagnostic=m_coupling->AuditSteadyState(candidate.current);
                m_coupling->WriteHistory(m_prefix+".steady-candidate.bin",candidate);
                const std::string path=m_prefix+".steady-candidate.json";Check(!std::filesystem::exists(path),"steady candidate evidence is immutable");std::ofstream out(path);Check(bool(out),"cannot create steady candidate evidence");
                out<<"{\"schema\":\"WJH27SteadyCandidateV1\",\"qualified_reference\":false,\"purpose\":\"replayable initial guess; no convergence or physical-time claim\",\"state\":";m_coupling->Report(out,candidate,diagnostic,true);out<<"}\n"<<std::flush;Check(bool(out),"steady candidate evidence write failed");out.close();Check(!out.fail(),"steady candidate evidence close failed");
            }
            Check(m_steadyResult.converged,"bounded native steady reference solve did not converge; original history preserved");
            m_coupling->Store(m_history.current);m_time=m_history.time;m_solveCompleted=true;return;
        }
        for(int n=0;n<m_steps;++n){const double start=Clock(),left=m_history.time,right=(m_history.step+1)*m_history.dt;
            WriteScientificFrame(TerminalReasons(left,right,true));
            m_last=m_coupling->Advance(m_history,int(m_history.step+1)==m_injectRetryStep);m_time=m_history.time;m_currentProcessAccepted=true;
            if(m_history.step%std::uint64_t(m_infoStride)==0||n+1==m_steps){m_coupling->Report(m_stepsLog,m_history,m_last,false);m_stepsLog<<'\n'<<std::flush;Check(bool(m_stepsLog),"coupled step-log write failed");std::cout<<"WJ-H27 step "<<m_history.step<<", t="<<m_time<<", coupled residual="<<m_last.residual.Maximum()<<", nonlinear iterations="<<m_last.nonlinearIterations<<std::endl;}
            if(m_recoveryStride>0&&m_history.step%std::uint64_t(m_recoveryStride)==0){double io=Clock();std::filesystem::path file=m_prefix+".recovery-"+std::to_string(m_history.step)+".bin";m_coupling->WriteHistory(file,m_history);m_recovery.push_back(file);
                while(m_recovery.size()>2){Check(std::filesystem::remove(m_recovery.front()),"owned coupled recovery rotation failed");m_removedRecovery.push_back(int(m_history.step));m_recovery.erase(m_recovery.begin());}m_recoverySeconds.push_back(Clock()-io);}
            unsigned reasons=TerminalReasons(left,right,false);if(m_scienceStride>0&&m_history.step%std::uint64_t(m_scienceStride)==0)reasons|=Periodic;if(m_terminalOutput&&n+1==m_steps)reasons|=Final;WriteScientificFrame(reasons);
            const double elapsed=Clock()-start;m_completeStepSeconds.push_back(elapsed);m_completeAdvanceSeconds+=elapsed;
        }
        m_coupling->Store(m_history.current);m_solveCompleted=true;
    }catch(const std::exception&e){
        try{
            m_coupling->WriteHistory(m_prefix+".failed-accepted.bin",m_history);
            WriteScientificIndex();
            auto diagnostic=m_history.step?m_last:m_coupling->AuditSteadyState(m_history.current,m_history.time);
            const std::string path=m_prefix+".failure.json";Check(!std::filesystem::exists(path),"coupled failure evidence is immutable");std::ofstream out(path);Check(bool(out),"cannot create coupled failure evidence");
            out<<std::setprecision(17)<<"{\"schema\":\"WJH27CoupledFailureV1\",\"status\":\"FAILED_WITH_ACCEPTED_HISTORY_PRESERVED\",\"reason\":"<<RANSWallSolveJson::Quote(e.what())<<",\"exception_type\":"<<RANSWallSolveJson::Quote(RANSWallSolveJson::ExceptionType(e))<<",\"elapsed_solve_seconds\":"<<Clock()-solveStart
               <<",\"accepted_state_residual_kind\":"<<RANSWallSolveJson::Quote(m_history.step?"last accepted complete time-step residual":"steady diagnostic of fresh initial state")<<",\"failed_attempt\":";AttemptJSON(out,m_coupling->LastAttemptInfo(),m_coupling->Options());out<<",\"state\":";m_coupling->Report(out,m_history,diagnostic,true);out<<"}\n"<<std::flush;Check(bool(out),"coupled failure evidence write failed");out.close();Check(!out.fail(),"coupled failure evidence close failed");
        }catch(const std::exception&journal){std::cerr<<"WJ-H27 failure-journal error: "<<journal.what()<<std::endl;}
        std::cerr<<"WJ-H27 coupled advancement: "<<e.what()<<std::endl;throw std::runtime_error(e.what());
    }
}
void WJNativeCoupled::v_Output()
{
    try{
        double start=Clock();if(m_terminalOutput&&m_steps==0)WriteScientificFrame(Final);
        m_coupling->WriteHistory(m_prefix+".final.bin",m_history);
        for(int i=0;i<5;++i){std::vector<Array<OneD,NekDouble>>coefficients{m_fields[i]->UpdateCoeffs()};std::vector<std::string>variables{m_session->GetVariable(i)};WriteFld(m_prefix+"-"+variables[0]+".fld",m_fields[i],coefficients,variables);}
        WriteScientificIndex();
        double exportSeconds=Clock()-start;std::ofstream report(m_prefix+".result.json");Check(bool(report),"cannot create coupled result");report<<std::setprecision(17)<<"{\"schema\":\"WJH27RegisteredCoupledRunV2\",\"factory\":\"WJNativeCoupledV2\",\"normal_main\":true,\"qualification\":\"execution evidence; closure gates are assessed independently\",\"construction_seconds\":"<<m_constructionSeconds<<",\"initialisation_seconds\":"<<m_initialisationSeconds<<",\"restart\":"<<(m_restart.empty()?"false":"true")<<",\"complete_advancing_seconds\":"<<m_completeAdvanceSeconds<<",\"final_export_seconds\":"<<exportSeconds<<",\"complete_step_seconds\":[";
        for(int i=0;i<int(m_completeStepSeconds.size());++i){if(i)report<<',';report<<m_completeStepSeconds[i];}report<<"],\"recovery_seconds\":[";for(int i=0;i<int(m_recoverySeconds.size());++i){if(i)report<<',';report<<m_recoverySeconds[i];}report<<"],\"rotated_recovery_count\":"<<m_removedRecovery.size()<<",\"retained_recovery\":[";for(int i=0;i<int(m_recovery.size());++i){if(i)report<<',';report<<'\"'<<m_recovery[i].string()<<'\"';}
        report << ']';
        if (!m_restart.empty())
        {
            report << ",\"recovery_load_seconds\":" << m_recoveryTiming.readSeconds
                   << ",\"recovery_cache_rebuild_seconds\":" << m_recoveryTiming.cacheRebuildSeconds;
        }
        report<<",\"scientific_index\":"<<RANSWallSolveJson::Quote(m_prefix+".science-index.json")<<",\"scientific_frame_seconds\":"<<m_scientificFrameSeconds<<",\"scientific_index_seconds\":"<<m_scientificIndexSeconds<<",\"scientific_timing_scope\":\"subsets of complete step/export totals; do not add twice\",\"steady_mode\":"<<(m_steadyMode?"true":"false")<<",\"steady_solver\":";if(m_steadyMode)SteadyJSON(report,m_steadyResult,m_steadyPreconditionerMass);else report<<"null";report<<",\"state\":";m_coupling->Report(report,m_history,m_last,m_coefficients);report<<"}\n"<<std::flush;Check(bool(report),"coupled result write failed");report.close();Check(!report.fail(),"coupled result close failed");
    }catch(const std::exception&e){std::cerr<<"WJ-H27 coupled output: "<<e.what()<<std::endl;throw std::runtime_error(e.what());}
}
void WJNativeCoupled::v_GenerateSummary(SolverUtils::SummaryList&summary)
{
    SolverUtils::AddSummaryItem(summary,"WJ coupled path","WJNativeCoupledV2; serial Q3/Q2, outer-only epsilon");
    SolverUtils::AddSummaryItem(summary,"Geometry",
        m_coupling->FlowKernel().Core().Mapping().normalMap ?
            "cosine hill normal layer and separate outer map" : "flat");
    SolverUtils::AddSummaryItem(summary,"Wall tensor",WJNative::WallModelV2);
    SolverUtils::AddSummaryItem(summary,"Integration",
        m_coupling->Options().UsesImplicitAdvection()
            ? "BDF2 with implicit advection, backward Euler startup"
            : "current coupled BDF2/EX2, BE startup, no clipping");
    SolverUtils::AddSummaryItem(summary,"Pressure gauge",m_coupling->FlowKernel().Gauge()?"native restricted mean":"fixed by total outlet traction");
}
void WJNativeCoupled::v_EvaluateExactSolution(unsigned int field,Array<OneD,NekDouble>&output,NekDouble time)
{
    const int n=m_fields[field]->GetTotPoints();output=Array<OneD,NekDouble>(n,0.);if(!m_coupling->Options().manufactured)return;
    const auto &exact=m_coupling->ManufacturedSolution();
    Array<OneD,NekDouble>x(n),z(n),d(n);m_fields[field]->GetCoords(x,z,d);
    for(int q=0;q<n;++q){if(field==3)output[q]=exact.OuterEpsilon(x[q],z[q],time).value;else{auto a=exact.Fields(x[q],z[q],time);output[q]=field==0?a.u.value:field==1?a.w.value:field==2?a.k.value:a.Pi.value;}}
}
}
