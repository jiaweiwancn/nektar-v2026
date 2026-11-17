// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeSteadyKrylov.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <Eigen/Core>
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
        {"projected_norm",J::Number(a.projectedNorm)},{"true_scaled_norm",J::Number(a.trueScaledNorm)},
        {"true_norm",J::Number(a.trueNorm)},{"true_checked",J::Boolean(a.trueChecked)},
        {"breakdown",J::Boolean(a.breakdown)},{"converged",J::Boolean(a.converged)}}));
    for(const auto&a:r.lineSearch)lineSearch.push_back(J::Object({{"newton",J::Number(a.newton)},
        {"trial",J::Number(a.trial)},{"step",J::Number(a.step)},{"true_norm",J::Number(a.trueNorm)},
        {"admissible",J::Boolean(a.admissible)},{"evaluated",J::Boolean(a.evaluated)},
        {"sufficient_decrease",J::Boolean(a.sufficientDecrease)},{"accepted",J::Boolean(a.accepted)},
        {"reason",J::Quote(a.reason)}}));
    for(const auto&a:r.newton)newton.push_back(J::Object({{"iteration",J::Number(a.iteration)},
        {"base_norm",J::Number(a.baseNorm)},{"forcing",J::Number(a.forcing)},
        {"linear_true_norm",J::Number(a.linearTrueNorm)},{"accepted_norm",J::Number(a.acceptedNorm)},
        {"accepted_step",J::Number(a.acceptedStep)},{"linear_converged",J::Boolean(a.linearConverged)},
        {"accepted",J::Boolean(a.accepted)}}));
    return J::Object({{"status",J::Quote(r.status)},{"reason",J::Quote(r.reason)},
        {"converged",J::Boolean(r.converged)},{"has_evaluated_state",J::Boolean(r.hasEvaluatedState)},
        {"best_from_accepted",J::Boolean(r.bestFromAcceptedIterate)},{"x",V(r.x)},{"residual",V(r.residual)},
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
std::string Probe(const std::string&name,const SteadyDerivativeProbe&p,const Vec&base,const Vec&direction,const Vec&exact)
{
    std::vector<std::string>actions;for(const auto&v:p.actions)actions.push_back(V(v));
    return J::Object({{"name",J::Quote(name)},{"completed",J::Boolean(p.completed)},{"reason",J::Quote(p.reason)},
        {"base",V(base)},{"direction",V(direction)},{"analytic_action",V(exact)},
        {"multipliers",J::Numbers(p.stepMultipliers)},{"actions",J::List(actions)},
        {"successive_scaled_changes",J::Numbers(p.successiveRelativeScaledChanges)},
        {"finite_differences",FD(p.finiteDifferences)},{"statistics",Statistics(p.statistics)}});
}
}
int main(int argc,char**argv)
{
    try
    {
        if(argc!=2)throw std::invalid_argument("usage: steady_krylov_probe new-output.json");
        if(std::filesystem::exists(argv[1]))throw std::invalid_argument("immutable probe output exists");
        std::vector<std::string>solves,probes,exceptions;
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
            SteadyKrylov solver(cb,Vector({1.}),Vector({2.}),options);const Vec initial=Vector({0.});
            const auto result=solver.Solve(initial);
            Check(!result.converged&&result.hasEvaluatedState&&!result.bestFromAcceptedIterate&&result.trueNorm<result.lastAcceptedNorm,
                "best unaccepted admissible state semantics were not exercised");
            solves.push_back(J::Object({{"name",J::Quote("best_unaccepted_state")},{"system",J::Quote("no_real_root")},
                {"initial",V(initial)},{"maximum_seconds",J::Number(options.maximumSeconds)},{"result",Result(result)}}));
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
        const auto result=J::Object({{"schema",J::Quote("WJH27SteadyKrylovProbeV1")},{"status",J::Quote("RECORDED")},
            {"solver_identity",J::Quote(SteadyKrylov::Identity)},{"unknown_scale",V(problem.sx)},
            {"physical_residual_scale",V(problem.b)},{"residual_row_scale",V(problem.RowScale())},
            {"acceptance_norm",J::Number(Options().acceptanceNorm)},
            {"scope",J::Quote("Standalone generic solver verification; no native coupled/PDE convergence or cost qualification")},
            {"solves",J::List(solves)},{"derivative_probes",J::List(probes)},{"fatal_exceptions",J::List(exceptions)}});
        std::ofstream output(argv[1]);output<<result<<'\n';Check(output.good(),"probe output write failed");
        std::cout<<"Recorded "<<solves.size()<<" solves, "<<probes.size()<<" derivative probes, "<<exceptions.size()<<" fatal controls\n";
    }
    catch(const std::exception&error){std::cerr<<error.what()<<'\n';return 2;}
    return 0;
}
