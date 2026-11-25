// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeSteadyKrylov.h>
#include <Eigen/QR>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <utility>

namespace Nektar::WJNative
{
namespace
{
using Clock=std::chrono::steady_clock;
using Vec=SteadyVector;
using Mat=Eigen::MatrixXd;

double Seconds(Clock::time_point start)
{return std::chrono::duration<double>(Clock::now()-start).count();}

struct Timer
{
    explicit Timer(double &destination):value(destination),start(Clock::now()){}
    ~Timer(){value+=Seconds(start);}
    double &value;
    Clock::time_point start;
};

void Require(bool condition,const std::string &message)
{if(!condition)throw std::invalid_argument(message);}

bool Positive(double x){return std::isfinite(x)&&x>0;}
bool Fraction(double x){return Positive(x)&&x<1;}

struct Stop : std::runtime_error
{explicit Stop(const std::string &message):std::runtime_error(message){}};

void CheckVector(const Vec &v,Eigen::Index size,const char *description)
{
    if(v.size()!=size||!v.allFinite())
        throw std::runtime_error(std::string(description)+" has incorrect size or nonfinite entries");
}

struct Evaluation
{
    bool admissible=false;
    bool nonlinearValidationAttempted=false, nonlinearValidated=false;
    Vec residual;
    std::string reason;
};

struct LinearResult
{
    Vec correction;
    bool converged=false;
    double trueNorm=0;
};

class Runner
{
public:
    Runner(const SteadyKrylovCallbacks &callbacks,const Vec &unknownScale,
           const Vec &rowScale,const SteadyKrylovOptions &options)
        :cb(callbacks),sx(unknownScale),sr(rowScale),o(options),start(Clock::now())
    {out.nonlinearValidationEnabled=bool(cb.ValidateNonlinearState);}

    SteadyKrylovResult Solve(const Vec &initial)
    {
        CheckVector(initial,sx.size(),"initial steady state");
        out.x=initial;
        try
        {
            Vec x=initial,f=Base(x);
            double norm=Norm(f),previousBase=norm;
            Observe(x,f,norm,true);Deadline();
            if(norm<=o.acceptanceNorm){Converged();return Finish();}
            for(int iteration=0;iteration<o.maximumNewtonIterations;++iteration)
            {
                newton=iteration;
                Deadline();
                // Re-evaluate accepted bases; no linearization/preconditioner
                // cache may substitute for the current nonlinear equations.
                if(iteration>0)
                {
                    f=Base(x);norm=Norm(f);Observe(x,f,norm,true);Deadline();
                    if(norm<=o.acceptanceNorm){Converged();return Finish();}
                }
                if(cb.PreparePreconditioner)
                {
                    ++out.statistics.preconditionerPreparations;
                    Timer timer(out.statistics.preparationSeconds);
                    cb.PreparePreconditioner(x,f);
                }
                Deadline();
                double forcing=o.forcingInitial;
                if(iteration>0)
                    forcing=std::clamp(o.forcingGamma*std::pow(norm/previousBase,o.forcingExponent),
                                       o.forcingMinimum,o.forcingMaximum);
                previousBase=norm;
                SteadyNewtonRecord record;
                record.iteration=iteration;record.baseNorm=norm;record.forcing=forcing;
                out.newton.push_back(record);
                const auto linear=Linear(x,f,norm,forcing);
                out.newton.back().linearTrueNorm=linear.trueNorm;
                out.newton.back().linearConverged=linear.converged;
                if(!linear.converged&&
                   (!o.allowIncompleteLinearStep||linear.trueNorm>=norm))
                    throw Stop("LINEAR_NOT_CONVERGED: declared Krylov cycles did not meet both current residual tests");
                const Vec correction=(sx.array()*linear.correction.array()).matrix();
                if(!correction.allFinite()||correction.cwiseAbs().maxCoeff()==0)
                    throw Stop("NEWTON_DIRECTION_UNUSABLE: zero or nonfinite correction");
                bool accepted=false;
                {
                    Timer timer(out.statistics.lineSearchSeconds);
                    double step=1;
                    for(int trial=0;trial<o.maximumLineSearchTrials;++trial,step*=o.lineSearchReduction)
                    {
                        Deadline();
                        SteadyLineSearchRecord ls;
                        ls.newton=iteration;ls.trial=trial;ls.step=step;
                        const Vec candidate=x+step*correction;
                        if(!candidate.allFinite())
                            throw Stop("NONFINITE_NEWTON_TRIAL: solver arithmetic failed; no clipping or hidden retry");
                        if((candidate.array()==x.array()).all())
                            throw Stop("NEWTON_ROUNDOFF: line-search correction no longer changes the state");
                        auto evaluation=Trial(candidate,false);
                        ls.admissible=evaluation.admissible;ls.reason=evaluation.reason;
                        ls.nonlinearValidationAttempted=evaluation.nonlinearValidationAttempted;
                        ls.nonlinearValidated=evaluation.nonlinearValidated;
                        if(evaluation.admissible)
                        {
                            ls.trueNorm=Norm(evaluation.residual);ls.evaluated=true;
                            Observe(candidate,evaluation.residual,ls.trueNorm,false);
                            // Norm-Armijo uses the caller's fixed physical norm.
                            ls.sufficientDecrease=ls.trueNorm<=o.acceptanceNorm||
                                ls.trueNorm<=(1-o.armijo*step)*norm;
                            if(!ls.sufficientDecrease)
                            {
                                ++out.statistics.lineSearchRejectedTrials;
                                ls.reason="current true norm fails sufficient decrease";
                            }
                        }
                        out.lineSearch.push_back(ls);Deadline();
                        if(ls.sufficientDecrease)
                        {
                            out.lineSearch.back().accepted=true;
                            x=candidate;f=std::move(evaluation.residual);norm=ls.trueNorm;
                            Observe(x,f,norm,true);accepted=true;
                            ++out.statistics.acceptedNewtonSteps;
                            auto &nr=out.newton.back();nr.accepted=true;
                            nr.acceptedNorm=norm;nr.acceptedStep=step;
                            break;
                        }
                    }
                }
                if(!accepted)throw Stop("LINE_SEARCH_NOT_CONVERGED: declared feasible Armijo trials exhausted");
                if(norm<=o.acceptanceNorm){Converged();return Finish();}
            }
            throw Stop("NEWTON_NOT_CONVERGED: declared Newton iteration bound exhausted");
        }
        catch(const Stop &error){out.reason=error.what();}
        return Finish();
    }

    SteadyDerivativeProbe Probe(const Vec &base,const Vec &direction,
                                const std::vector<double> &multipliers)
    {
        CheckVector(base,sx.size(),"derivative probe base");
        CheckVector(direction,sx.size(),"derivative probe direction");
        Require(!multipliers.empty(),"derivative probe needs prescribed step multipliers");
        for(double multiplier:multipliers)
            Require(Positive(multiplier)&&Positive(o.fdRelativeStep*multiplier)&&
                    o.fdRelativeStep*multiplier>=o.fdMinimumRelativeStep,
                    "invalid derivative probe step multiplier");
        SteadyDerivativeProbe result;result.stepMultipliers=multipliers;
        result.nonlinearValidationEnabled=bool(cb.ValidateNonlinearState);
        try
        {
            const Vec f=Base(base);
            const Vec scaledDirection=(direction.array()/sx.array()).matrix();
            for(double multiplier:multipliers)
            {
                Vec action=Jacobian(base,f,scaledDirection,multiplier,"direction-probe",-1,-1);
                if(!result.actions.empty())
                {
                    const Vec a=(sr.array()*action.array()).matrix();
                    const Vec b=(sr.array()*result.actions.back().array()).matrix();
                    const double denominator=std::max(a.stableNorm(),b.stableNorm());
                    result.successiveRelativeScaledChanges.push_back(
                        denominator==0?0:(a-b).stableNorm()/denominator);
                }
                result.actions.push_back(std::move(action));Deadline();
            }
            result.completed=true;result.reason="fixed prescribed directional probes completed; accuracy not certified";
        }
        catch(const Stop &error){result.reason=error.what();}
        out.statistics.totalSeconds=Seconds(start);
        result.statistics=out.statistics;
        result.finiteDifferences=std::move(out.finiteDifferences);
        return result;
    }

private:
    void Deadline()const
    {if(Seconds(start)>=o.maximumSeconds)throw Stop("DEADLINE: finite steady-solve wall budget exhausted");}

    SteadyTrialValidity Valid(const Vec &x)
    {
        Deadline();++out.statistics.validityEvaluations;
        Timer timer(out.statistics.validitySeconds);
        auto value=cb.ValidTrial(x);
        if(!value.admissible&&value.reason.empty())
            throw std::runtime_error("ValidTrial rejected a state without an explanatory reason");
        return value;
    }

    Vec Residual(const Vec &x)
    {
        Deadline();++out.statistics.residualEvaluations;
        Timer timer(out.statistics.residualSeconds);
        Vec f=cb.Residual(x);CheckVector(f,sx.size(),"steady residual callback");return f;
    }

    void ValidateNonlinear(const Vec &x)
    {
        if(!cb.ValidateNonlinearState)return;
        Deadline();++out.statistics.nonlinearValidationEvaluations;
        Timer timer(out.statistics.nonlinearValidationSeconds);
        cb.ValidateNonlinearState(x);
    }

    Vec Base(const Vec &x)
    {
        try
        {
            const auto validity=Valid(x);
            if(!validity.admissible)
                throw std::runtime_error("inadmissible steady base: "+validity.reason);
            ValidateNonlinear(x);
            return Residual(x);
        }
        catch(const SteadyInadmissibleTrial &error)
        {throw std::runtime_error(std::string("inadmissible steady base: ")+error.what());}
    }

    Evaluation Trial(const Vec &x,bool finiteDifference)
    {
        Evaluation result;
        try
        {
            const auto validity=Valid(x);
            result.admissible=validity.admissible;result.reason=validity.reason;
            if(result.admissible)
            {
                if(!finiteDifference&&cb.ValidateNonlinearState)
                {
                    result.nonlinearValidationAttempted=true;
                    ValidateNonlinear(x);
                    result.nonlinearValidated=true;
                }
                result.residual=Residual(x);
            }
        }
        catch(const SteadyInadmissibleTrial &error)
        {
            result.admissible=false;result.reason=error.what();
            if(result.reason.empty())throw std::runtime_error("empty inadmissible-trial reason");
            if(result.nonlinearValidationAttempted&&!result.nonlinearValidated)
                ++out.statistics.nonlinearValidationRejectedTrials;
        }
        if(!result.admissible)
        {
            if(finiteDifference)++out.statistics.fdRejectedTrials;
            else ++out.statistics.lineSearchRejectedTrials;
        }
        return result;
    }

    double Norm(const Vec &f)
    {
        Deadline();++out.statistics.normEvaluations;
        Timer timer(out.statistics.normSeconds);
        const double value=cb.TrueNorm(f);
        if(!std::isfinite(value)||value<0)
            throw std::runtime_error("TrueNorm returned a negative or nonfinite value");
        return value;
    }

    void Observe(const Vec &x,const Vec &f,double norm,bool accepted)
    {
        if(!out.hasEvaluatedState||norm<out.trueNorm)
        {
            out.hasEvaluatedState=true;out.x=x;out.residual=f;out.trueNorm=norm;out.bestNewtonIteration=newton;
            out.bestFromAcceptedIterate=accepted;
        }
        else if(accepted&&(x.array()==out.x.array()).all())out.bestFromAcceptedIterate=true;
        if(accepted)out.lastAcceptedNorm=norm;
    }

    void Converged()
    {out.converged=true;out.status="CONVERGED";out.reason="current true steady residual meets caller acceptance norm";}

    SteadyKrylovResult Finish()
    {out.statistics.totalSeconds=Seconds(start);return std::move(out);}

    Vec Precondition(const Vec &scaledResidual)
    {
        Deadline();++out.statistics.preconditionerApplications;
        Timer timer(out.statistics.preconditionerSeconds);
        const Vec rhs=(scaledResidual.array()/sr.array()).matrix();
        CheckVector(rhs,sx.size(),"unscaled preconditioner RHS");
        const Vec physical=cb.Precondition(rhs);
        CheckVector(physical,sx.size(),"preconditioner callback correction");
        Vec scaled=(physical.array()/sx.array()).matrix();
        CheckVector(scaled,sx.size(),"scaled preconditioner correction");return scaled;
    }

    // Returns raw physical J(x)*diag(sx)*scaledDirection. Both sides must be
    // explicitly tested before choosing a central or feasible one-sided rule.
    Vec Jacobian(const Vec &x,const Vec &f,const Vec &scaledDirection,
                 double multiplier,const std::string &purpose,int cycle,int iteration)
    {
        Deadline();++out.statistics.jacobianActions;
        CheckVector(scaledDirection,sx.size(),"scaled Jacobian direction");
        const double directionNorm=scaledDirection.cwiseAbs().maxCoeff();
        if(directionNorm==0)return Vec::Zero(sx.size());
        const Vec y=(x.array()/sx.array()).matrix();
        const Vec direction=(sx.array()*scaledDirection.array()).matrix();
        CheckVector(y,sx.size(),"scaled steady base");
        CheckVector(direction,sx.size(),"physical Jacobian direction");
        const double coordinateNorm=std::max(1.,y.cwiseAbs().maxCoeff());
        double relative=o.fdRelativeStep*multiplier;
        for(int reduction=0;reduction<=o.maximumFDReductions;++reduction)
        {
            Deadline();
            const auto attemptStart=Clock::now();
            const double displacement=relative*coordinateNorm,h=displacement/directionNorm;
            if(!Positive(h)||!Positive(displacement))
                throw Stop("FD_ARITHMETIC: nonfinite or vanishing prescribed displacement");
            const Vec plus=x+h*direction,minus=x-h*direction;
            if(!plus.allFinite()||!minus.allFinite())
                throw Stop("FD_ARITHMETIC: nonfinite trial; no unreported adaptive step");
            if((plus.array()==x.array()).all()||(minus.array()==x.array()).all())
                throw Stop("FD_ROUNDOFF: prescribed displacement is not representable");
            SteadyFDRecord record;record.newton=newton;record.cycle=cycle;record.iteration=iteration;
            record.reduction=reduction;record.purpose=purpose;record.relativeStep=relative;
            record.parameterStep=h;record.scaledCoordinateStep=displacement;
            record.scaledDirectionInfinityNorm=directionNorm;
            out.finiteDifferences.push_back(record);
            auto &audit=out.finiteDifferences.back();
            auto positive=Trial(plus,true);
            audit.plusAdmissible=positive.admissible;audit.plusReason=positive.reason;
            auto negative=Trial(minus,true);
            audit.minusAdmissible=negative.admissible;audit.minusReason=negative.reason;
            Vec action;
            if(positive.admissible&&negative.admissible)
            {
                action=(positive.residual-negative.residual)/(2*h);audit.method="central";
            }
            else if(positive.admissible||negative.admissible)
            {
                const bool forward=positive.admissible;
                const Vec half=x+(forward?.5:-.5)*h*direction;
                if((half.array()==x.array()).all())
                    throw Stop("FD_ROUNDOFF: feasible one-sided half step is not representable");
                auto middle=Trial(half,true);
                audit.halfAdmissible=middle.admissible;audit.halfReason=middle.reason;
                if(middle.admissible)
                {
                    // Difference from base first to avoid a large -3*f term.
                    if(forward)action=(4*(middle.residual-f)-(positive.residual-f))/h;
                    else action=(4*(f-middle.residual)-(f-negative.residual))/h;
                    audit.method=forward?"forward-second-order":"backward-second-order";
                }
            }
            audit.seconds=Seconds(attemptStart);
            if(action.size()!=0)
            {
                if(cb.CorrectionAction)
                {
                    const Vec correction=cb.CorrectionAction(x,direction);
                    CheckVector(correction,sx.size(),"additive correction action");
                    action+=correction;
                }
                if(!action.allFinite())throw Stop("FD_ARITHMETIC: nonfinite directional derivative");
                Deadline();return action;
            }
            audit.method="infeasible-reduction";
            // This path is reachable only after actual callback infeasibility.
            if(reduction==o.maximumFDReductions||relative*o.fdReduction<o.fdMinimumRelativeStep)
                throw Stop("FD_INFEASIBLE: declared derivative reductions cannot produce an admissible stencil");
            relative*=o.fdReduction;
        }
        throw Stop("FD_INFEASIBLE: derivative construction exhausted");
    }

    LinearResult Linear(const Vec &x,const Vec &f,double baseNorm,double forcing)
    {
        Timer timer(out.statistics.linearSeconds);
        const Vec b=-(sr.array()*f.array()).matrix();
        CheckVector(b,sx.size(),"scaled Newton RHS");
        const double bNorm=b.stableNorm();
        if(!Positive(bNorm))throw Stop("LINEAR_SCALING: nonzero true residual has zero/nonfinite scaled norm");
        const int n=static_cast<int>(sx.size()),m=std::min(n,o.restartDimension);
        Vec best=Vec::Zero(n);double bestNorm=baseNorm,bestScaled=bNorm;
        for(int cycle=0;cycle<o.maximumLinearCycles;++cycle)
        {
            Deadline();++out.statistics.linearCycles;
            const Vec baseCorrection=best;
            Vec r=b;
            if(baseCorrection.cwiseAbs().maxCoeff()!=0)
                r-=(sr.array()*Jacobian(x,f,baseCorrection,1.,"restart-residual",cycle,-1).array()).matrix();
            CheckVector(r,n,"current scaled linear residual");
            const double beta=r.stableNorm();
            const double initialTrue=Norm((r.array()/sr.array()).matrix());
            SteadyLinearRecord restart;restart.newton=newton;restart.cycle=cycle;
            restart.iteration=-1;restart.trueChecked=true;restart.projectedNorm=beta;
            restart.trueScaledNorm=beta;restart.trueNorm=initialTrue;
            restart.converged=beta<=forcing*bNorm&&initialTrue<=forcing*baseNorm;
            out.linear.push_back(restart);
            if(restart.converged)
                return {baseCorrection,true,initialTrue};
            if(!Positive(beta))throw Stop("LINEAR_SCALING: inconsistent zero/nonfinite Krylov residual");
            Mat v=Mat::Zero(n,m+1),z=Mat::Zero(n,m),h=Mat::Zero(m+1,m);
            v.col(0)=r/beta;
            for(int j=0;j<m;++j)
            {
                Deadline();++out.statistics.linearIterations;
                z.col(j)=Precondition(v.col(j));
                Vec w=(sr.array()*Jacobian(x,f,z.col(j),1.,"arnoldi",cycle,j).array()).matrix();
                CheckVector(w,n,"scaled Arnoldi action");
                const double before=w.stableNorm();
                if(!std::isfinite(before))throw Stop("ARNOLDI_ARITHMETIC: nonfinite action norm");
                double after=0;
                {
                    Timer orthogonalization(out.statistics.arnoldiSeconds);
                    // Reorthogonalized modified Gram-Schmidt, including the
                    // second-pass corrections in the original Hessenberg.
                    for(int pass=0;pass<2;++pass)
                        for(int i=0;i<=j;++i)
                        {
                            const double coefficient=v.col(i).dot(w);
                            h(i,j)+=coefficient;w-=coefficient*v.col(i);
                        }
                    after=w.stableNorm();h(j+1,j)=after;
                }
                if(!std::isfinite(after)||!h.col(j).allFinite())
                    throw Stop("ARNOLDI_ARITHMETIC: nonfinite reorthogonalization");
                const bool breakdown=after<=o.arnoldiBreakdownRelative*before;
                if(!breakdown)v.col(j+1)=w/after;
                const int k=j+1;
                const Mat small=h.topLeftCorner(k+1,k);
                Vec target=Vec::Zero(k+1);target[0]=beta;
                // Rank-revealing small dense QR also handles happy/singular
                // Arnoldi breakdown without dividing by a zero diagonal.
                Eigen::CompleteOrthogonalDecomposition<Mat> qr(small);
                const Vec coefficients=qr.solve(target);
                CheckVector(coefficients,k,"small Krylov least-squares solution");
                SteadyLinearRecord record;record.newton=newton;record.cycle=cycle;
                record.iteration=j;record.leastSquaresRank=static_cast<int>(qr.rank());
                record.projectedNorm=(target-small*coefficients).stableNorm();record.breakdown=breakdown;
                if(breakdown||k==m||k%o.trueLinearCheckInterval==0||record.projectedNorm<=forcing*bNorm)
                {
                    const Vec candidate=baseCorrection+z.leftCols(k)*coefficients;
                    // This is a new nonlinear-residual FD evaluation in the
                    // accumulated direction, not the Arnoldi recurrence norm.
                    const Vec check=b-(sr.array()*Jacobian(x,f,candidate,1.,"true-linear-residual",cycle,j).array()).matrix();
                    CheckVector(check,n,"recomputed linear residual");
                    record.trueChecked=true;record.trueScaledNorm=check.stableNorm();
                    record.trueNorm=Norm((check.array()/sr.array()).matrix());
                    record.converged=record.trueScaledNorm<=forcing*bNorm&&record.trueNorm<=forcing*baseNorm;
                    if(record.trueNorm<bestNorm||
                       (record.trueNorm==bestNorm&&record.trueScaledNorm<bestScaled))
                    {best=candidate;bestNorm=record.trueNorm;bestScaled=record.trueScaledNorm;}
                    out.linear.push_back(record);Deadline();
                    if(record.converged)return {candidate,true,record.trueNorm};
                }
                else out.linear.push_back(record);
                if(breakdown)break;
            }
        }
        return {best,false,bestNorm};
    }

    const SteadyKrylovCallbacks &cb;
    const Vec &sx,&sr;
    const SteadyKrylovOptions &o;
    Clock::time_point start;
    int newton=-1;
    SteadyKrylovResult out;
};
}

SteadyKrylov::SteadyKrylov(SteadyKrylovCallbacks callbacks,SteadyVector unknownScale,
    SteadyVector residualRowScale,SteadyKrylovOptions options)
    :m_callbacks(std::move(callbacks)),m_unknownScale(std::move(unknownScale)),
     m_residualRowScale(std::move(residualRowScale)),m_options(options)
{
    const auto &o=m_options;
    Require(bool(m_callbacks.Residual)&&bool(m_callbacks.TrueNorm)&&
            bool(m_callbacks.Precondition)&&bool(m_callbacks.ValidTrial),
            "steady solver requires residual, fixed true norm, physical preconditioner and feasibility callbacks");
    Require(m_unknownScale.size()>0&&m_unknownScale.size()<std::numeric_limits<int>::max()&&
            m_unknownScale.size()==m_residualRowScale.size()&&m_unknownScale.allFinite()&&
            m_residualRowScale.allFinite()&&(m_unknownScale.array()>0).all()&&
            (m_residualRowScale.array()>0).all(),"steady solver scales must be matching positive finite vectors");
    Require(Positive(o.maximumSeconds)&&Positive(o.acceptanceNorm),
            "steady solver needs explicit finite wall budget and true-norm acceptance threshold");
    Require(o.maximumNewtonIterations>0&&o.restartDimension>0&&o.maximumLinearCycles>0&&
            o.trueLinearCheckInterval>0&&o.maximumLineSearchTrials>0&&o.maximumFDReductions>=0,
            "invalid declared steady algorithm cycle bounds");
    Require(Positive(o.fdRelativeStep)&&Positive(o.fdMinimumRelativeStep)&&
            o.fdMinimumRelativeStep<=o.fdRelativeStep&&Fraction(o.fdReduction)&&
            Fraction(o.lineSearchReduction)&&Fraction(o.armijo)&&o.armijo<.5,
            "invalid declared FD or norm-Armijo policy");
    Require(Fraction(o.forcingInitial)&&Fraction(o.forcingMinimum)&&Fraction(o.forcingMaximum)&&
            o.forcingMinimum<=o.forcingInitial&&o.forcingInitial<=o.forcingMaximum&&
            Positive(o.forcingGamma)&&Positive(o.forcingExponent)&&Fraction(o.arnoldiBreakdownRelative),
            "invalid declared inexact-Newton forcing or Arnoldi policy");
}

SteadyKrylovResult SteadyKrylov::Solve(const SteadyVector &initial)const
{return Runner(m_callbacks,m_unknownScale,m_residualRowScale,m_options).Solve(initial);}

SteadyDerivativeProbe SteadyKrylov::ProbeDirection(const SteadyVector &base,
    const SteadyVector &direction,const std::vector<double> &multipliers)const
{return Runner(m_callbacks,m_unknownScale,m_residualRowScale,m_options).Probe(base,direction,multipliers);}
}
