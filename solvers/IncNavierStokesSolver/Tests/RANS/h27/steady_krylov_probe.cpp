// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeSteadyKrylov.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <Eigen/Core>
#include <Eigen/LU>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <thread>

using namespace Nektar::WJNative;
namespace J=Nektar::RANSWallSolveJson;
using Vec=SteadyVector;
namespace
{
void Check(bool value,const std::string&message){if(!value)throw std::runtime_error(message);}
Vec Vector(std::initializer_list<double>values)
{Vec v(values.size());int i=0;for(double x:values)v[i++]=x;return v;}
std::string V(const Vec&v)
{std::vector<double>values;values.reserve(v.size());for(Eigen::Index i=0;i<v.size();++i)values.push_back(v[i]);return J::Numbers(values);}
std::string Statistics(const SteadyKrylovStatistics&s)
{
    return J::Object({{"residual_evaluations",J::Number(s.residualEvaluations)},
        {"validity_evaluations",J::Number(s.validityEvaluations)},{"norm_evaluations",J::Number(s.normEvaluations)},
        {"preconditioner_applications",J::Number(s.preconditionerApplications)},
        {"preconditioner_preparations",J::Number(s.preconditionerPreparations)},
        {"jacobian_actions",J::Number(s.jacobianActions)},{"fd_rejections",J::Number(s.fdRejectedTrials)},
        {"line_search_rejections",J::Number(s.lineSearchRejectedTrials)},
        {"linear_iterations",J::Number(s.linearIterations)},{"linear_cycles",J::Number(s.linearCycles)},
        {"accepted_newton_steps",J::Number(s.acceptedNewtonSteps)},
        {"nonlinear_validation_evaluations",J::Number(s.nonlinearValidationEvaluations)},
        {"nonlinear_validation_rejections",J::Number(s.nonlinearValidationRejectedTrials)},
        {"nonlinear_validation_seconds",J::Number(s.nonlinearValidationSeconds)},
        {"residual_seconds",J::Number(s.residualSeconds)},{"validity_seconds",J::Number(s.validitySeconds)},
        {"norm_seconds",J::Number(s.normSeconds)},{"preconditioner_seconds",J::Number(s.preconditionerSeconds)},
        {"preparation_seconds",J::Number(s.preparationSeconds)},{"arnoldi_seconds",J::Number(s.arnoldiSeconds)},
        {"linear_seconds",J::Number(s.linearSeconds)},{"line_search_seconds",J::Number(s.lineSearchSeconds)},
        {"total_seconds",J::Number(s.totalSeconds)}});
}
std::string FD(const std::vector<SteadyFDRecord>&records)
{
    std::vector<std::string>out;for(const auto&r:records)out.push_back(J::Object({
        {"newton",J::Number(r.newton)},{"cycle",J::Number(r.cycle)},{"iteration",J::Number(r.iteration)},
        {"reduction",J::Number(r.reduction)},{"purpose",J::Quote(r.purpose)},{"method",J::Quote(r.method)},
        {"relative_step",J::Number(r.relativeStep)},{"parameter_step",J::Number(r.parameterStep)},
        {"scaled_coordinate_step",J::Number(r.scaledCoordinateStep)},
        {"scaled_direction_norm",J::Number(r.scaledDirectionInfinityNorm)},
        {"plus_admissible",J::Boolean(r.plusAdmissible)},{"minus_admissible",J::Boolean(r.minusAdmissible)},
        {"half_admissible",J::Boolean(r.halfAdmissible)},{"plus_reason",J::Quote(r.plusReason)},
        {"minus_reason",J::Quote(r.minusReason)},{"half_reason",J::Quote(r.halfReason)}}));return J::List(out);
}
std::string Result(const SteadyKrylovResult&r)
{
    std::vector<std::string>linear,lineSearch,newton;
    for(const auto&a:r.linear)linear.push_back(J::Object({{"newton",J::Number(a.newton)},
        {"cycle",J::Number(a.cycle)},{"iteration",J::Number(a.iteration)},{"rank",J::Number(a.leastSquaresRank)},
        {"projected_norm",J::Number(a.projectedNorm)},
        {"true_scaled_norm",a.trueChecked?J::Number(a.trueScaledNorm):"null"},
        {"true_norm",a.trueChecked?J::Number(a.trueNorm):"null"},{"true_checked",J::Boolean(a.trueChecked)},
        {"breakdown",J::Boolean(a.breakdown)},{"converged",J::Boolean(a.converged)}}));
    for(const auto&a:r.lineSearch)lineSearch.push_back(J::Object({{"newton",J::Number(a.newton)},
        {"trial",J::Number(a.trial)},{"step",J::Number(a.step)},{"true_norm",J::Number(a.trueNorm)},
        {"admissible",J::Boolean(a.admissible)},{"evaluated",J::Boolean(a.evaluated)},
        {"predicted_norm",a.predictedNorm?J::Number(*a.predictedNorm):"null"},
        {"actual_to_predicted_decrease",a.decreaseRatio?J::Number(*a.decreaseRatio):"null"},
        {"nonlinear_validation_attempted",J::Boolean(a.nonlinearValidationAttempted)},
        {"nonlinear_validated",J::Boolean(a.nonlinearValidated)},
        {"sufficient_decrease",J::Boolean(a.sufficientDecrease)},{"accepted",J::Boolean(a.accepted)},
        {"reason",J::Quote(a.reason)}}));
    for(const auto&a:r.newton)newton.push_back(J::Object({{"iteration",J::Number(a.iteration)},
        {"base_norm",J::Number(a.baseNorm)},{"forcing",J::Number(a.forcing)},
        {"linear_true_norm",a.linearEvaluated?J::Number(a.linearTrueNorm):"null"},
        {"accepted_norm",J::Number(a.acceptedNorm)},
        {"accepted_step",J::Number(a.acceptedStep)},{"linear_converged",J::Boolean(a.linearConverged)},
        {"accepted",J::Boolean(a.accepted)}}));
    return J::Object({{"status",J::Quote(r.status)},{"reason",J::Quote(r.reason)},
        {"converged",J::Boolean(r.converged)},{"has_evaluated_state",J::Boolean(r.hasEvaluatedState)},
        {"nonlinear_validation_enabled",J::Boolean(r.nonlinearValidationEnabled)},
        {"best_from_accepted",J::Boolean(r.bestFromAcceptedIterate)},{"x",V(r.x)},{"residual",V(r.residual)},
        {"last_accepted_x",V(r.lastAcceptedX)},
        {"true_norm",J::Number(r.trueNorm)},{"last_accepted_norm",J::Number(r.lastAcceptedNorm)},
        {"statistics",Statistics(r.statistics)},{"finite_differences",FD(r.finiteDifferences)},
        {"linear",J::List(linear)},{"line_search",J::List(lineSearch)},{"newton",J::List(newton)}});
}
struct Problem
{
    Vec sx=Vector({1e-3,3.,2e3,.02}),b=Vector({1e4,.002,7.,.3});
    Vec rowWeight=Vector({.25,2.,.5,1.5}),prepared=Vec::Ones(4);
    Eigen::Matrix4d M;
    std::vector<std::string>preconditionerSamples;
    Problem(){M<<1,.3,-.2,.1,.2,1,.4,-.1,-.1,.3,1,.2,.4,-.2,.1,1;}
    Vec X(const Vec&y)const{return (sx.array()*y.array()).matrix();}
    Vec Y(const Vec&x)const{return (x.array()/sx.array()).matrix();}
    Vec RowScale()const{return (rowWeight.array()/b.array()).matrix();}
    Vec F(const Vec&x)const
    {
        const Vec y=Y(x);
        if(y[1]<=0)throw SteadyInadmissibleTrial("typed residual log-domain rejection");
        const Vec f=Vector({y[0]*y[0]-1,std::log(y[1]),y[2]*y[2]*y[2]-1,std::exp(y[3])-std::exp(1.)});
        return (b.array()*(M*f).array()).matrix();
    }
    Vec Diagonal(const Vec&y)const
    {return Vector({2*y[0],1/y[1],3*y[2]*y[2],std::exp(y[3])});}
    Vec Action(const Vec&x,const Vec&direction)const
    {const Vec product=(Diagonal(Y(x)).array()*Y(direction).array()).matrix();return (b.array()*(M*product).array()).matrix();}
    double Norm(const Vec&f)const{return (f.array()/b.array()).abs().maxCoeff();}
    SteadyKrylovCallbacks Callbacks(bool typed=false)
    {
        SteadyKrylovCallbacks cb;
        cb.Residual=[this](const Vec&x){return F(x);};
        cb.TrueNorm=[this](const Vec&f){return Norm(f);};
        cb.ValidTrial=[this,typed](const Vec&x)
        {
            const auto y=Y(x);
            if(!typed&&(y[0]<=0||y[1]<=0||y[2]<=0))return SteadyTrialValidity{false,"positive nonlinear coordinate required"};
            return SteadyTrialValidity{};
        };
        cb.PreparePreconditioner=[this](const Vec&x,const Vec&f)
        {Check((F(x)-f).norm()==0,"prepared base residual is stale");prepared=Y(x);};
        cb.Precondition=[this](const Vec&rhs)
        {
            const Vec correction=(sx.array()*rhs.array()/b.array()/Diagonal(prepared).array()).matrix();
            if(preconditionerSamples.size()<12)preconditionerSamples.push_back(J::Object({
                {"prepared_y",V(prepared)},{"raw_rhs",V(rhs)},{"physical_correction",V(correction)}}));
            return correction;
        };
        return cb;
    }
};
SteadyKrylovOptions Options()
{
    SteadyKrylovOptions o;o.maximumSeconds=5;o.acceptanceNorm=1e-10;
    o.restartDimension=4;o.maximumNewtonIterations=40;o.maximumLinearCycles=3;
    o.trueLinearCheckInterval=1;o.forcingInitial=o.forcingMinimum=o.forcingMaximum=.001;
    o.fdRelativeStep=1e-5;return o;
}

std::vector<std::string> OriginalMeritControls()
{
    std::vector<std::string> receipts;
    auto solve=[](const Eigen::Matrix2d &a,const Vec &rhs,const Vec &x,
                  SteadyKrylovOptions o,bool approximate,int &observed) {
        double mass=0.;
        Eigen::Matrix2d inverse=Eigen::Matrix2d::Identity();
        SteadyKrylovCallbacks cb;
        cb.Residual=[&](const Vec &v)->Vec{return rhs+a*v;};
        cb.ValidTrial=[](const Vec&){return SteadyTrialValidity{};};
        cb.TrueNorm=[](const Vec &v){return v.lpNorm<Eigen::Infinity>();};
        cb.ResidualBlocks=[](const Vec &v){
            return std::array<double,6>{std::abs(v[0]),std::abs(v[1]),0.,0.,0.,0.};
        };
        cb.SetCorrectionMass=[&](double value){mass=value;};
        cb.CorrectionAction=[&](const Vec&,const Vec &d)->Vec{return mass*d;};
        cb.PreparePreconditioner=[&](const Vec&,const Vec&) {
            if(approximate)inverse<<2.,0.,-.4,1.;
            else inverse=(a+mass*Eigen::Matrix2d::Identity()).inverse();
        };
        cb.Precondition=[&](const Vec &v)->Vec{return inverse*v;};
        cb.ObserveAccepted=[&](const Vec&,double){++observed;};
        return SteadyKrylov(cb,Vec::Ones(2),Vec::Ones(2),o).Solve(x);
    };
    auto record=[&](const char *name,const SteadyKrylovResult &r) {
        std::vector<std::string> directions;
        for(const auto &d:r.meritDirections)
            directions.push_back(J::Object({
                {"mass",J::Number(d.mass)},{"forcing",J::Number(d.forcing)},
                {"error",J::Number(d.linearNorm)},{"scaled_error",J::Number(d.scaledError)},
                {"spacing_uncertainty",J::Number(d.actionUncertainty)},
                {"slope",J::Number(d.activeSlope)},{"step",J::Number(d.step)},
                {"blocks",J::Numbers(d.activeBlocks)},{"decision",J::Quote(d.decision)}}));
        receipts.push_back(J::Object({{"name",J::Quote(name)},{"result",Result(r)},
            {"directions",J::List(directions)},{"controller_start",r.meritStart->ToJson()},
            {"controller_end",r.meritEnd->ToJson()}}));
    };
    auto o=Options();
    o.maximumNewtonIterations=1;
    o.forcingInitial=.5;
    o.merit.emplace();
    o.merit->maximumWork=1;
    o.merit->initialMass=3.;
    const Vec zero=Vec::Zero(2);
    Eigen::Matrix2d a;
    a<<-2.,0.,0.,1.;
    int seen=0;
    const auto ascent=solve(a,Vector({1,0}),zero,o,false,seen);
    Check(ascent.meritDirections.size()==1&&ascent.meritDirections[0].decision=="RESOLVED_ASCENT"&&
          ascent.meritDirections[0].linearNorm<1e-8&&ascent.lineSearch.empty()&&seen==1&&
          ascent.lastAcceptedX==zero,"exact shifted ascent reached an actual trial");
    record("exact_shifted_ascent",ascent);

    o.merit->initialMass=1.;
    a<<-2.,0.,0.,-.5;
    seen=0;
    const auto tie=solve(a,Vector({1,1}),zero,o,false,seen);
    Check(tie.meritDirections[0].decision=="RESOLVED_ASCENT"&&
          tie.meritDirections[0].activeBlocks.size()==2&&tie.lineSearch.empty()&&seen==1,
          "tied original merit used a favorable block");
    record("active_block_tie",tie);

    a<<0.,.005,1.,-1.;
    seen=0;
    const auto small=solve(a,Vector({1,0}),zero,o,false,seen);
    Check(small.meritDirections[0].decision=="WEAK_VERIFIED_DESCENT"&&
          small.meritDirections[0].step<1./64.&&small.newton[0].accepted&&
          small.lastAcceptedNorm<1.&&seen==2,"small verified descent was lost across active-block change");
    record("finite_secant_caveat",small);

    o.merit->maximumWork=12;
    a<<.1,3.,0.,.2;
    seen=0;
    const auto tighter=solve(a,Vector({1,0}),zero,o,true,seen);
    Check(tighter.meritDirections.size()>=2&&tighter.meritDirections[0].decision=="RESOLVED_ASCENT"&&
          tighter.meritDirections[0].forcing==.5&&tighter.meritDirections[1].forcing==.1&&
          tighter.meritDirections[0].mass==tighter.meritDirections[1].mass&&
          tighter.inverses[0].serial==tighter.inverses[1].serial&&
          tighter.meritDirections[1].linearNorm<tighter.meritDirections[0].linearNorm&&
          tighter.newton[0].accepted&&seen==2,"same-base same-inverse tightening control failed");
    record("same_base_tighter",tighter);

    a<<-2.,0.,0.,1.;
    o.merit->initialMass=3.;
    seen=0;
    const auto changed=solve(a,Vector({1,0}),zero,o,false,seen);
    Check(changed.meritDirections.size()>1&&changed.meritDirections[0].decision=="RESOLVED_ASCENT"&&
          changed.meritDirections.back().mass!=3.&&changed.newton[0].accepted&&seen==2&&
          changed.lastAcceptedNorm<1e-8,"changed mass did not check a fresh direction at the same base");
    record("changed_mass",changed);

    a<<4.,0.,0.,1.;
    o.merit->initialMass=1.;
    seen=0;
    const auto interrupted=solve(a,Vector({1,0}),zero,o,false,seen);
    const auto saved=interrupted.meritEnd;
    Check(saved->goodSteps==1&&saved->nextMass==1.,"controller did not preserve first useful mass");
    auto continuous=o;
    continuous.maximumNewtonIterations=2;
    const auto uninterrupted=solve(a,Vector({1,0}),zero,continuous,false,seen);
    o.merit->restored=saved;
    const auto resumed=solve(a,Vector({1,0}),interrupted.lastAcceptedX,o,false,seen);
    Check(resumed.meritStart->ToJson()==saved->ToJson()&&
          resumed.meritEnd->ToJson()==uninterrupted.meritEnd->ToJson()&&
          resumed.meritEnd->nextMass==.5&&resumed.lastAcceptedX==uninterrupted.lastAcceptedX,
          "controller resume changed effective mass or next decision for identical action");
    record("controller_resume",resumed);
    return receipts;
}

std::string Probe(const std::string&name,const SteadyDerivativeProbe&p,const Vec&base,const Vec&direction,const Vec&exact)
{
    std::vector<std::string>actions;for(const auto&v:p.actions)actions.push_back(V(v));
    return J::Object({{"name",J::Quote(name)},{"completed",J::Boolean(p.completed)},{"reason",J::Quote(p.reason)},
        {"nonlinear_validation_enabled",J::Boolean(p.nonlinearValidationEnabled)},
        {"base",V(base)},{"direction",V(direction)},{"analytic_action",V(exact)},
        {"multipliers",J::Numbers(p.stepMultipliers)},{"actions",J::List(actions)},
        {"successive_scaled_changes",J::Numbers(p.successiveRelativeScaledChanges)},
        {"finite_differences",FD(p.finiteDifferences)},{"statistics",Statistics(p.statistics)}});
}
SteadyKrylovCallbacks ScalarCallbacks()
{
    SteadyKrylovCallbacks cb;
    cb.Residual=[](const Vec&x){return Vector({x[0]-1});};
    cb.TrueNorm=[](const Vec&f){return std::abs(f[0]);};
    cb.Precondition=[](const Vec&rhs){return Vec(rhs);};
    cb.ValidTrial=[](const Vec&x){return SteadyTrialValidity{std::abs(x[0])<=2,"scalar constitutive domain"};};
    return cb;
}
}
int main(int argc,char**argv)
{
    try
    {
        if(argc!=2)throw std::invalid_argument("usage: steady_krylov_probe new-output.json");
        if(std::filesystem::exists(argv[1]))throw std::invalid_argument("immutable probe output exists");
        std::vector<std::string>solves,probes,exceptions,validationControls;
        const Vec initialY=Vector({.3,4.,.4,1.8});
        for(bool typed:{false,true})
        {
            Problem problem;auto options=Options();const auto initial=problem.X(initialY);const Vec copy=initial;
            SteadyKrylov solver(problem.Callbacks(typed),problem.sx,problem.RowScale(),options);
            const auto result=solver.Solve(initial);
            Check((copy.array()==initial.array()).all(),"Solve mutated caller's initial vector");
            Check(result.converged&&result.trueNorm<=options.acceptanceNorm,"analytic coupled solve did not converge");
            Check(result.statistics.lineSearchRejectedTrials>0,"feasible Newton line search was not exercised");
            Check(result.statistics.preconditionerPreparations>0&&result.statistics.preconditionerApplications>0,"right preconditioner was not exercised");
            solves.push_back(J::Object({{"name",J::Quote(typed?"typed_residual_rejection":"validity_rejection")},
                {"system",J::Quote("coupled_analytic")},{"maximum_seconds",J::Number(options.maximumSeconds)},
                {"initial",V(initial)},{"result",Result(result)},
                {"preconditioner_samples",J::List(problem.preconditionerSamples)}}));
        }
        for(const std::string name:{"central","forward","backward","typed_forward","reduced_central","infeasible_stencil"})
        {
            Problem problem;auto cb=problem.Callbacks();auto options=Options();
            const Vec y=Vector({.8,name=="central"?1.3:1.,.7,1.1});const Vec base=problem.X(y);
            const Vec direction=problem.X(Vector({.3,.4,.2,.1}));
            if(name=="forward"||name=="backward"||name=="reduced_central"||name=="infeasible_stencil")
                cb.ValidTrial=[&problem,name](const Vec&x)
                {
                    const double z=problem.Y(x)[1];bool ok=name=="forward"?z>=1:name=="backward"?z<=1:
                        name=="reduced_central"?std::abs(z-1)<=1.5e-6:z==1;
                    return SteadyTrialValidity{ok,ok?"":"explicit derivative-domain boundary"};
                };
            if(name=="typed_forward")
            {
                cb.ValidTrial=[](const Vec&){return SteadyTrialValidity{};};
                cb.Residual=[&problem](const Vec&x)
                {if(problem.Y(x)[1]<1)throw SteadyInadmissibleTrial("typed lower-bound stencil rejection");return problem.F(x);};
            }
            if(name=="infeasible_stencil")options.maximumFDReductions=2;
            SteadyKrylov solver(cb,problem.sx,problem.RowScale(),options);
            const auto result=solver.ProbeDirection(base,direction);
            Check(result.completed==(name!="infeasible_stencil"),"derivative completion status differs");
            probes.push_back(Probe(name,result,base,direction,problem.Action(base,direction)));
        }
        {
            auto cb=ScalarCallbacks();
            auto options=Options();
            options.maximumNewtonIterations=1;
            options.measureDirection=true;
            int additiveCalls=0;
            double frozenBase=std::numeric_limits<double>::quiet_NaN();
            cb.PrepareCorrection=[&](const Vec&x,const Vec&){frozenBase=x[0];};
            cb.CorrectionAction=[&](const Vec&x,const Vec&direction){
                Check(x[0]==frozenBase,"additive correction applied at a different nonlinear base");
                ++additiveCalls;
                return Vector({2.*direction[0]});
            };
            cb.Precondition=[](const Vec&rhs){return Vec(rhs/3.);};
            const Vec one=Vector({1.});
            SteadyKrylov solver(cb,one,one,options);
            const auto fresh=solver.ProbeDirection(Vector({.125}),one);
            Check(fresh.completed&&frozenBase==.125&&std::abs(fresh.actions.back()[0]-3.)<1e-9,
                  "standalone shifted probe did not prepare its first base");
            const auto result=solver.Solve(Vector({0.}));
            Check(!result.converged&&result.hasEvaluatedState&&std::abs(result.x[0]-1./3)<1e-10&&
                  std::abs(result.trueNorm-2./3)<1e-10,"correction action changed the nonlinear residual");
            Check(!result.newton.empty()&&result.newton[0].linearConverged&&additiveCalls>=2,
                  "Arnoldi or true linear check skipped the correction action");
            const auto &trial=result.lineSearch.at(0);
            Check(trial.predictedNorm&&std::abs(*trial.predictedNorm-2./3)<1e-10&&
                  trial.decreaseRatio&&std::abs(*trial.decreaseRatio-1.)<1e-9,
                  "shifted correction did not predict the original residual");
            const auto probe=solver.ProbeDirection(Vector({.25}),one);
            Check(probe.completed&&frozenBase==.25&&std::abs(probe.actions.back()[0]-3.)<1e-9,
                  "additive derivative action disagrees with the analytic shifted derivative");
            solves.push_back(J::Object({{"name",J::Quote("additive_correction_only")},
                {"system",J::Quote("scalar_analytic")},{"result",Result(result)}}));
            probes.push_back(Probe("additive_correction_only",probe,Vector({.25}),one,Vector({3.})));
        }
        for(const std::string name:{"newton_bound","linear_bound","line_search_bound","deadline_before_base","deadline_after_base"})
        {
            Problem problem;auto cb=problem.Callbacks();auto options=Options();
            if(name=="newton_bound")options.maximumNewtonIterations=1;
            if(name=="linear_bound"){options.restartDimension=1;options.maximumLinearCycles=1;}
            if(name=="line_search_bound")options.maximumLineSearchTrials=1;
            if(name=="deadline_before_base")options.maximumSeconds=std::numeric_limits<double>::min();
            if(name=="deadline_after_base")
            {
                options.maximumSeconds=.05;auto prepare=cb.PreparePreconditioner;
                cb.PreparePreconditioner=[prepare](const Vec&x,const Vec&f)
                {prepare(x,f);std::this_thread::sleep_for(std::chrono::milliseconds(75));};
            }
            SteadyKrylov solver(cb,problem.sx,problem.RowScale(),options);const auto initial=problem.X(initialY);
            const auto result=solver.Solve(initial);Check(!result.converged&&result.status=="NOT_CONVERGED","bounded control falsely converged");
            solves.push_back(J::Object({{"name",J::Quote(name)},{"system",J::Quote("coupled_analytic")},
                {"initial",V(initial)},{"maximum_seconds",J::Number(options.maximumSeconds)},{"result",Result(result)}}));
        }
        {
            // No real root. The full Newton trial improves the residual but
            // fails deliberately strict Armijo, exercising returned BEST vs
            // accepted-state semantics without mutating the caller's state.
            SteadyKrylovCallbacks cb;cb.Residual=[](const Vec&x){return Vector({1-x[0]+.99*x[0]*x[0]});};
            cb.TrueNorm=[](const Vec&f){return std::abs(f[0]);};
            cb.Precondition=[](const Vec&rhs){return Vec(-rhs);};cb.ValidTrial=[](const Vec&){return SteadyTrialValidity{};};
            auto options=Options();options.armijo=.4;options.maximumLineSearchTrials=1;
            std::vector<double> observed;
            cb.ObserveAccepted=[&](const Vec &x,double){observed.push_back(x[0]);};
            SteadyKrylov solver(cb,Vector({1.}),Vector({2.}),options);const Vec initial=Vector({0.});
            const auto result=solver.Solve(initial);
            Check(!result.converged&&result.hasEvaluatedState&&!result.bestFromAcceptedIterate&&result.trueNorm<result.lastAcceptedNorm,
                "best unaccepted admissible state semantics were not exercised");
            Check(observed==std::vector<double>{initial[0]},
                  "the uncommitted trial was published to the accepted observer");
            solves.push_back(J::Object({{"name",J::Quote("best_unaccepted_state")},{"system",J::Quote("no_real_root")},
                {"initial",V(initial)},{"maximum_seconds",J::Number(options.maximumSeconds)},{"result",Result(result)}}));
        }
        for(const std::string name:{"rejected_candidate_not_best","guarded_linesearch_progress"})
        {
            auto cb=ScalarCallbacks();auto options=Options();options.maximumNewtonIterations=1;
            options.maximumLineSearchTrials=name=="rejected_candidate_not_best"?1:2;
            int forbiddenResiduals=0;std::vector<double>validated;
            cb.Residual=[&](const Vec&x){if(x[0]>.75)++forbiddenResiduals;return Vector({x[0]-1});};
            cb.ValidateNonlinearState=[&](const Vec&x)
            {validated.push_back(x[0]);if(x[0]>.75)throw SteadyInadmissibleTrial("full-principal-candidate-marker");};
            const auto result=SteadyKrylov(cb,Vector({1}),Vector({2}),options).Solve(Vector({0}));
            Check(!result.converged&&result.hasEvaluatedState&&result.bestFromAcceptedIterate,
                  "full-guard rejection contaminated returned-state provenance");
            Check(forbiddenResiduals==0&&result.statistics.nonlinearValidationRejectedTrials==1,
                  "full-guard-rejected candidate was residual-evaluated or not counted");
            Check(std::abs(result.x[0]-(name=="rejected_candidate_not_best"?0.:.5))<1e-8,
                  "full-guard rejected point became best or accepted state");
            validationControls.push_back(J::Object({{"name",J::Quote(name)},{"result",Result(result)},
                {"full_validation_inputs",J::Numbers(validated)},
                {"forbidden_candidate_residual_evaluations",J::Number(forbiddenResiduals)}}));
        }
        for(const std::string name:{"fd_outside_full_margin","fd_domain_guard","fd_typed_domain_guard"})
        {
            auto cb=ScalarCallbacks();std::vector<double>validated,residualInputs;
            cb.ValidTrial=[&](const Vec&x){const bool ok=name!="fd_domain_guard"||x[0]>=0;
                return SteadyTrialValidity{ok,ok?"":"FD-domain-boundary-marker"};};
            cb.Residual=[&](const Vec&x){residualInputs.push_back(x[0]);
                if(name=="fd_typed_domain_guard"&&x[0]<0)throw SteadyInadmissibleTrial("typed-FD-domain-marker");
                return Vector({x[0]+x[0]*x[0]*x[0]});};
            cb.ValidateNonlinearState=[&](const Vec&x){validated.push_back(x[0]);
                if(std::abs(x[0])>1e-12)throw SteadyInadmissibleTrial("nonlinear-full-margin-marker");};
            const auto result=SteadyKrylov(cb,Vector({1}),Vector({2}),Options()).ProbeDirection(Vector({0}),Vector({1}));
            Check(result.completed&&validated.size()==1&&validated[0]==0,
                  "full nonlinear guard leaked into FD probes or did not guard the probe base");
            Check(result.statistics.nonlinearValidationEvaluations==1&&result.statistics.nonlinearValidationRejectedTrials==0,
                  "wrong full-validation provenance on FD probe");
            validationControls.push_back(J::Object({{"name",J::Quote(name)},
                {"probe",Probe(name,result,Vector({0}),Vector({1}),Vector({1}))},
                {"full_validation_inputs",J::Numbers(validated)},{"residual_inputs",J::Numbers(residualInputs)}}));
        }
        {
            auto cb=ScalarCallbacks();std::vector<double>validated;
            cb.Residual=[](const Vec&x){return Vector({1-x[0]+.99*x[0]*x[0]});};
            cb.Precondition=[](const Vec&rhs){return Vec(-rhs);};
            cb.ValidateNonlinearState=[&](const Vec&x){validated.push_back(x[0]);};
            auto options=Options();options.armijo=.4;options.maximumLineSearchTrials=1;
            const auto result=SteadyKrylov(cb,Vector({1}),Vector({2}),options).Solve(Vector({0}));
            Check(!result.converged&&result.hasEvaluatedState&&!result.bestFromAcceptedIterate&&
                  result.trueNorm<result.lastAcceptedNorm&&validated.size()==2,
                  "fully validated best unaccepted candidate was lost");
            validationControls.push_back(J::Object({{"name",J::Quote("guarded_best_unaccepted_state")},
                {"result",Result(result)},{"full_validation_inputs",J::Numbers(validated)}}));
        }
        {
            Problem baseline,guarded;auto options=Options();std::vector<std::string>validated;
            auto cb=guarded.Callbacks();cb.ValidateNonlinearState=[&](const Vec&x){validated.push_back(V(x));};
            const auto a=SteadyKrylov(baseline.Callbacks(),baseline.sx,baseline.RowScale(),options).Solve(baseline.X(initialY));
            const auto b=SteadyKrylov(cb,guarded.sx,guarded.RowScale(),options).Solve(guarded.X(initialY));
            Check(a.converged&&b.converged&&(a.x.array()==b.x.array()).all()&&(a.residual.array()==b.residual.array()).all(),
                  "no-op full validation changed the established nonlinear algorithm");
            validationControls.push_back(J::Object({{"name",J::Quote("noop_full_guard_equivalence")},
                {"baseline",Result(a)},{"guarded",Result(b)},{"full_validation_inputs",J::List(validated)}}));
        }
        auto mustThrow=[&](const std::string&name,const std::string&marker,const auto&operation)
        {
            try{operation();}catch(const std::exception&error)
            {Check(std::string(error.what()).find(marker)!=std::string::npos,"wrong fatal error for "+name+": "+error.what());
             exceptions.push_back(J::Object({{"name",J::Quote(name)},{"reason",J::Quote(error.what())}}));return;}
            throw std::runtime_error("fatal callback/input error was swallowed: "+name);
        };
        Problem problem;const Vec base=problem.X(Vector({.8,1.,.7,1.1})),direction=problem.X(Vector({.3,.4,.2,.1}));
        mustThrow("fatal_trial_residual","fatal-residual-marker",[&]
        {
            auto cb=problem.Callbacks();cb.Residual=[&](const Vec&x)
            {if((x.array()!=base.array()).any())throw std::logic_error("fatal-residual-marker");return problem.F(x);};
            SteadyKrylov(cb,problem.sx,problem.RowScale(),Options()).ProbeDirection(base,direction);
        });
        mustThrow("fatal_trial_validity","fatal-validity-marker",[&]
        {
            auto cb=problem.Callbacks();cb.ValidTrial=[&](const Vec&x)
            {if((x.array()!=base.array()).any())throw std::logic_error("fatal-validity-marker");return SteadyTrialValidity{};};
            SteadyKrylov(cb,problem.sx,problem.RowScale(),Options()).ProbeDirection(base,direction);
        });
        mustThrow("fatal_preconditioner","fatal-preconditioner-marker",[&]
        {
            auto cb=problem.Callbacks();cb.Precondition=[](const Vec&)->Vec{throw std::logic_error("fatal-preconditioner-marker");};
            SteadyKrylov(cb,problem.sx,problem.RowScale(),Options()).Solve(problem.X(initialY));
        });
        mustThrow("typed_base_is_fatal","inadmissible steady base",[&]
        {
            auto cb=problem.Callbacks();cb.Residual=[](const Vec&)->Vec{throw SteadyInadmissibleTrial("base-domain-marker");};
            SteadyKrylov(cb,problem.sx,problem.RowScale(),Options()).Solve(base);
        });
        mustThrow("malformed_residual","incorrect size",[&]
        {auto cb=problem.Callbacks();cb.Residual=[](const Vec&){return Vec::Zero(3);};SteadyKrylov(cb,problem.sx,problem.RowScale(),Options()).Solve(base);});
        mustThrow("invalid_true_norm","TrueNorm returned",[&]
        {auto cb=problem.Callbacks();cb.TrueNorm=[](const Vec&){return -1.;};SteadyKrylov(cb,problem.sx,problem.RowScale(),Options()).Solve(base);});
        mustThrow("missing_wall_cap","explicit finite wall budget",[&]
        {auto options=Options();options.maximumSeconds=0;SteadyKrylov(problem.Callbacks(),problem.sx,problem.RowScale(),options);});
        mustThrow("invalid_unknown_scale","positive finite vectors",[&]
        {Vec bad=problem.sx;bad[0]=-1;SteadyKrylov(problem.Callbacks(),bad,problem.RowScale(),Options());});
        int initialGuardResiduals=0;
        mustThrow("typed_full_initial_is_fatal","inadmissible steady base: full-initial-marker",[&]
        {
            auto cb=ScalarCallbacks();cb.Residual=[&](const Vec&x){++initialGuardResiduals;return Vector({x[0]-1});};
            cb.ValidateNonlinearState=[](const Vec&){throw SteadyInadmissibleTrial("full-initial-marker");};
            SteadyKrylov(cb,Vector({1}),Vector({2}),Options()).Solve(Vector({0}));
        });
        Check(initialGuardResiduals==0,"failed full initial guard reached Residual/Observe");
        mustThrow("typed_full_probe_base_is_fatal","inadmissible steady base: full-probe-base-marker",[&]
        {
            auto cb=ScalarCallbacks();cb.ValidateNonlinearState=[](const Vec&){throw SteadyInadmissibleTrial("full-probe-base-marker");};
            SteadyKrylov(cb,Vector({1}),Vector({2}),Options()).ProbeDirection(Vector({0}),Vector({1}));
        });
        mustThrow("fatal_full_candidate","fatal-full-candidate-marker",[&]
        {
            auto cb=ScalarCallbacks();cb.ValidateNonlinearState=[](const Vec&x){
                if(x[0]>.75)throw std::logic_error("fatal-full-candidate-marker");};
            SteadyKrylov(cb,Vector({1}),Vector({2}),Options()).Solve(Vector({0}));
        });
        int fullBaseCalls=0;
        mustThrow("typed_full_revalidated_base_is_fatal","inadmissible steady base: full-revalidated-base-marker",[&]
        {
            auto cb=ScalarCallbacks();double previous=std::numeric_limits<double>::quiet_NaN();
            cb.Residual=[](const Vec&x){return Vector({std::exp(x[0])-2});};
            cb.ValidateNonlinearState=[&](const Vec&x){++fullBaseCalls;
                if(x[0]==previous)throw SteadyInadmissibleTrial("full-revalidated-base-marker");previous=x[0];};
            SteadyKrylov(cb,Vector({1}),Vector({2}),Options()).Solve(Vector({0}));
        });
        Check(fullBaseCalls==3,"accepted nonlinear base was not fully revalidated without a generic cache");
        {
            auto cb=ScalarCallbacks();auto options=Options();
            cb.Residual=[](const Vec &x){return Vector({x[0]*x[0]*x[0]-1.});};
            cb.TrueNorm=[](const Vec &f) {
                if(f[0]>.3&&f[0]<.5)
                    std::this_thread::sleep_for(std::chrono::milliseconds(120));
                return std::abs(f[0]);
            };
            options.armijo=.4;options.maximumSeconds=.08;
            const Vec start=Vector({.72}),one=Vector({1.});
            const auto stopped=SteadyKrylov(cb,one,one,options).Solve(start);
            Check(stopped.reason.rfind("DEADLINE",0)==0&&!stopped.bestFromAcceptedIterate&&
                  stopped.lastAcceptedX.size()==1&&stopped.lastAcceptedX[0]==start[0]&&stopped.x[0]!=start[0]&&
                  stopped.trueNorm<stopped.lastAcceptedNorm,
                  "deadline replaced the last accepted state with a better uncommitted trial");
            options.maximumSeconds=3.;
            cb.TrueNorm=[](const Vec &f){return std::abs(f[0]);};
            const auto resumed=SteadyKrylov(cb,one,one,options).Solve(stopped.lastAcceptedX);
            Check(resumed.converged&&std::abs(resumed.x[0]-1.)<1e-10,"accepted-state resume failed");
            solves.push_back(J::Object({{"name",J::Quote("deadline_uncommitted_then_resume")},
                {"result",Result(stopped)},{"resumed",Result(resumed)}}));
        }
        {
            Problem p;auto options=Options();
            options.maximumInverseAge=4;options.reuseLinearIterations=1;
            const auto solver=SteadyKrylov(p.Callbacks(),p.sx,p.RowScale(),options);
            const auto result=solver.Solve(p.X(initialY));
            const auto again=solver.Solve(p.X(initialY));
            int refreshes=0;
            for(const auto &use:result.inverses)
                if(use.rebuild&&use.reason=="current linear check or work cap")++refreshes;
            Check(refreshes>0&&result.converged&&again.converged&&again.inverses.at(0).rebuild&&
                  again.inverses.at(0).age==0,
                  "reused inverse work cap or cold entry was not enforced");
            solves.push_back(J::Object({{"name",J::Quote("reuse_forced_cold_refresh")},
                {"refreshes",J::Number(refreshes)},{"result",Result(result)}}));
        }
        const auto result=J::Object({{"schema",J::Quote("WJH27SteadyKrylovProbeV1")},{"status",J::Quote("RECORDED")},
            {"solver_identity",J::Quote(SteadyKrylov::Identity)},{"unknown_scale",V(problem.sx)},
            {"physical_residual_scale",V(problem.b)},{"residual_row_scale",V(problem.RowScale())},
            {"acceptance_norm",J::Number(Options().acceptanceNorm)},
            {"scope",J::Quote("Standalone generic solver verification; no native coupled/PDE convergence or cost qualification")},
            {"solves",J::List(solves)},{"derivative_probes",J::List(probes)},{"fatal_exceptions",J::List(exceptions)},
            {"nonlinear_validation_controls",J::List(validationControls)},
            {"original_merit_controls",J::List(OriginalMeritControls())},
            {"full_initial_guard_residual_calls",J::Number(initialGuardResiduals)},
            {"revalidated_full_base_calls",J::Number(fullBaseCalls)}});
        std::ofstream output(argv[1]);output<<result<<'\n';Check(output.good(),"probe output write failed");
        std::cout<<"Recorded "<<solves.size()<<" solves, "<<probes.size()<<" derivative probes, "<<exceptions.size()<<" fatal controls\n";
    }
    catch(const std::exception&error){std::cerr<<error.what()<<'\n';return 2;}
    return 0;
}
