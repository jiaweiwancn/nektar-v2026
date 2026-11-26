// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeSteadyKrylov.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
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
    Vec action;
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
            if(o.merit)StartMerit(norm);
            Observe(x,f,norm,true);Deadline();
            if(norm<=o.acceptanceNorm){Converged();return Finish();}
            if(o.evaluateInitialOnly)
            {
                out.status="EVALUATED";
                out.reason="initial residual evaluated without applying a correction";
                return Finish();
            }
            int factorBase=-1;
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
                double forcing=o.forcingInitial;
                if(iteration>0)
                    forcing=std::clamp(o.forcingGamma*std::pow(norm/previousBase,o.forcingExponent),
                                       o.forcingMinimum,o.forcingMaximum);
                if(o.merit)forcing=out.meritEnd->nextForcing;
                previousBase=norm;
                SteadyNewtonRecord record;
                record.iteration=iteration;record.baseNorm=norm;record.forcing=forcing;
                record.missingLinearReason="correction base preparation";
                out.newton.push_back(record);
                if(cb.PrepareCorrection)cb.PrepareCorrection(x,f);
                if(o.merit)
                    RecoverDirection(x,f,norm,factorBase);
                else
                {
                    const bool reuse=factorBase>=0&&iteration-factorBase<=o.maximumInverseAge;
                    auto linear=Correct(x,f,norm,forcing,factorBase,
                        reuse?"reuse":factorBase<0?"cold entry":"age cap",reuse);
                    if(reuse&&!linear.converged)
                        linear=Correct(x,f,norm,forcing,factorBase,"current linear check or work cap",false);
                    RecordLinear(linear,forcing);
                    if(!linear.converged&&
                       (!o.allowIncompleteLinearStep||linear.trueNorm>=norm))
                        throw Stop("LINEAR_NOT_CONVERGED: declared Krylov cycles did not meet both current residual tests");
                    const Vec correction=PhysicalCorrection(linear);
                    Vec originalAction;
                    double modelSlope=0.;
                    if(o.measureDirection)
                    {
                        originalAction=std::move(linear.action);
                        if(cb.CorrectionAction)originalAction-=cb.CorrectionAction(x,correction);
                        if(o.adaptiveTrialStep)
                        {
                            constexpr double localFraction=.015625;
                            modelSlope=(Norm(f+localFraction*originalAction)-norm)/localFraction;
                        }
                    }
                    if(!Backtrack(x,f,norm,correction,originalAction,1.,modelSlope,o.maximumLineSearchTrials))
                        throw Stop("LINE_SEARCH_NOT_CONVERGED: declared feasible Armijo trials exhausted");
                }
                if(norm<=o.acceptanceNorm){Converged();return Finish();}
            }
            throw Stop("NEWTON_NOT_CONVERGED: declared Newton iteration bound exhausted");
        }
        catch(const Stop &error){out.reason=error.what();}
        catch(const SteadyCorrectionFailure &error){out.reason=error.what();}
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
            if(cb.PrepareCorrection)cb.PrepareCorrection(base,f);
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
        catch(const SteadyCorrectionFailure &error){result.reason=error.what();}
        out.statistics.totalSeconds=Seconds(start);
        result.statistics=out.statistics;
        result.finiteDifferences=std::move(out.finiteDifferences);
        return result;
    }

private:

    void StartMerit(double norm)
    {
        SteadyMeritState state;
        if(o.merit->restored)
        {
            state=*o.merit->restored;
            if(std::abs(state.acceptedNorm-norm)>2e-12*std::max(norm,state.acceptedNorm)+1e-15)
                throw std::invalid_argument("restored original merit differs from the accepted packet");
        }
        else
        {
            state.usedMass=state.nextMass=o.merit->initialMass;
            state.referenceNorm=state.previousNorm=state.acceptedNorm=norm;
            state.usedForcing=state.nextForcing=o.forcingInitial;
            state.decision="initial";
        }
        out.meritStart=out.meritEnd=state;
        if(cb.ObserveMeritState)cb.ObserveMeritState(state);
    }

    double StrongerMass(double mass)const
    {
        const double cap=std::max({1.,o.merit->initialMass,o.merit->reanchorMass});
        return mass>cap/2.?cap:2.*mass;
    }

    void AdvanceMerit(double baseNorm,const SteadyLineSearchRecord &trial)
    {
        auto state=*out.meritEnd;
        state.goodSteps=state.usedMass==currentMass?state.goodSteps:0;
        state.usedMass=state.nextMass=currentMass;
        state.usedForcing=out.newton.back().forcing;
        state.previousNorm=baseNorm;
        state.acceptedNorm=trial.trueNorm;
        state.step=trial.step;
        state.ratio=trial.decreaseRatio;
        ++state.accepted;
        const double relativeDecrease=(baseNorm-trial.trueNorm)/baseNorm;
        const bool reliable=trial.step>=.75&&trial.decreaseRatio&&
                            *trial.decreaseRatio>=.7&&*trial.decreaseRatio<=1.5&&relativeDecrease>=.2;
        state.goodSteps=reliable?state.goodSteps+1:0;
        state.decision="hold";
        if(state.goodSteps==2&&currentMass>0.)
        {
            state.nextMass=.5*currentMass;
            state.goodSteps=0;
            state.decision="weaken after two well-predicted large reductions";
        }
        else if(trial.step<.25&&relativeDecrease<.03&&StrongerMass(currentMass)>currentMass)
        {
            state.nextMass=StrongerMass(currentMass);
            state.decision="strengthen after limited actual progress";
        }
        state.nextForcing=std::clamp(o.forcingGamma*std::pow(trial.trueNorm/baseNorm,o.forcingExponent),
                                     o.forcingMinimum,o.forcingMaximum);
        out.meritEnd=state;
        out.meritAcceptances.push_back(state);
        if(cb.ObserveMeritState)cb.ObserveMeritState(state);
    }

    void RecordLinear(const LinearResult &linear,double forcing)
    {
        auto &entry=out.newton.back();
        entry.forcing=forcing;
        entry.linearEvaluated=true;
        entry.missingLinearReason.clear();
        entry.linearTrueNorm=linear.trueNorm;
        entry.linearConverged=linear.converged;
    }

    Vec PhysicalCorrection(const LinearResult &linear)const
    {
        Vec dx=sx.cwiseProduct(linear.correction);
        if(!dx.allFinite()||dx.cwiseAbs().maxCoeff()==0.)
            throw Stop("NEWTON_DIRECTION_UNUSABLE: zero or nonfinite correction");
        return dx;
    }

    // Evaluate the original model using the current checked linear action.
    // A second, predetermined spacing supplies a sensitivity allowance.
    // Reuse the primary action for the decision regardless of the second sign.
    SteadyMeritRecord MeasureMerit(const Vec &x,const Vec &f,double norm,
                                  const LinearResult &linear,double forcing,
                                  int proposal,int usedWork,Vec &original)
    {
        SteadyMeritRecord rec;
        rec.base=newton;
        rec.proposal=proposal;
        rec.inverse=int(out.inverses.size())-1;
        rec.work=usedWork;
        rec.mass=currentMass;
        rec.forcing=forcing;
        rec.baseNorm=norm;
        rec.baseBlocks=cb.ResidualBlocks(f);
        const bool nonzero=linear.correction.cwiseAbs().maxCoeff()!=0.;
        const Vec primary=nonzero?Jacobian(x,f,linear.correction,1.,"original-merit-primary",-1,proposal):
                                  linear.action;
        const double primarySpacing=nonzero?out.finiteDifferences.back().parameterStep:0.;
        const Vec error=f+primary;
        rec.linearNorm=Norm(error);
        rec.scaledError=sr.cwiseProduct(error).stableNorm();
        rec.scaledBaseNorm=sr.cwiseProduct(f).stableNorm();
        rec.shiftedErrorBlocks=cb.ResidualBlocks(error);
        rec.linearPassed=linear.converged&&rec.linearNorm<=forcing*norm&&
                         rec.scaledError<=forcing*rec.scaledBaseNorm;
        if(!nonzero)
        {
            rec.decision="LINEAR_UNRESOLVED";
            return rec;
        }
        const Vec shift=cb.CorrectionAction(x,PhysicalCorrection(linear));
        original=primary-shift;
        rec.correctionBlocks=cb.ResidualBlocks(shift);
        rec.fullPredictedBlocks=cb.ResidualBlocks(error-shift);
        const double multiplier=std::max(.5,o.fdMinimumRelativeStep/o.fdRelativeStep);
        rec.decision="ACTION_UNRESOLVED";
        if(multiplier>=1.)return rec;
        const Vec finer=Jacobian(x,f,linear.correction,multiplier,"original-merit-spacing",-1,proposal);
        if(out.finiteDifferences.back().parameterStep>=.8*primarySpacing)return rec;
        const Vec change=finer-primary;
        rec.actionUncertainty=2.*Norm(change);
        rec.scaledUncertainty=2.*sr.cwiseProduct(change).stableNorm();
        rec.modelEvaluated=true;
        const double actionNorm=Norm(original);
        const double roundoff=64.*std::numeric_limits<double>::epsilon()*std::max(norm,actionNorm);
        rec.slopeUncertainty=rec.actionUncertainty+roundoff;
        if(std::abs(*std::max_element(rec.baseBlocks.begin(),rec.baseBlocks.end())-norm)>roundoff)
            throw std::invalid_argument("original merit is not the maximum of its declared blocks");
        // Polarization gives each active Hilbert-seminorm slope without relying
        // on a finite one-sided secant across a change in the active block.
        rec.activeSlope=-std::numeric_limits<double>::infinity();
        const double t=actionNorm>norm?norm/actionNorm:1.;
        const auto plus=cb.ResidualBlocks(f+t*original);
        const auto minus=cb.ResidualBlocks(f-t*original);
        for(std::size_t i=0;i<rec.baseBlocks.size();++i)
            if(rec.baseBlocks[i]>=norm-64.*std::numeric_limits<double>::epsilon()*norm)
            {
                rec.activeBlocks.push_back(int(i));
                if(rec.baseBlocks[i]==norm)
                {
                    const double slope=(plus[i]-minus[i])/(4.*t)*
                        ((plus[i]+minus[i])/rec.baseBlocks[i]);
                    rec.activeSlope=std::max(rec.activeSlope,slope);
                }
            }
        if(!std::isfinite(rec.activeSlope))
            throw std::invalid_argument("original merit has no evaluable active block slope");

        const double gamma=std::max(.01,2.*o.armijo);
        for(double step=1.;step>=1./4096.;step*=.5)
        {
            const Vec prediction=f+step*original;
            const double value=Norm(prediction);
            if(value+step*rec.actionUncertainty+roundoff<=(1.-gamma*step)*norm)
            {
                rec.step=step;
                rec.predictedNorm=value;
                rec.predictedBlocks=cb.ResidualBlocks(prediction);
                break;
            }
        }
        if(!rec.linearPassed)
            rec.decision="LINEAR_UNRESOLVED";
        else if(rec.actionUncertainty>=.5*forcing*norm||
                rec.scaledUncertainty>=.5*forcing*rec.scaledBaseNorm)
            rec.decision="ACTION_UNRESOLVED";
        else if(rec.step>0.)
            rec.decision=(norm-rec.predictedNorm)>=.01*norm?"USEFUL_DESCENT":"WEAK_VERIFIED_DESCENT";
        else if(rec.activeSlope>rec.slopeUncertainty)
            rec.decision="RESOLVED_ASCENT";
        else if(rec.activeSlope<-rec.slopeUncertainty)
            rec.decision="INSUFFICIENT_MODEL_DECREASE";
        return rec;
    }

    void RecoverDirection(Vec &x,Vec &f,double &norm,int &factorBase)
    {
        const double initialMass=out.meritEnd->nextMass;
        currentMass=initialMass;
        bool refreshed=false;
        const auto initialWork=out.statistics.linearIterations;
        auto remaining=[&] {return o.merit->maximumWork-int(out.statistics.linearIterations-initialWork);};
        std::string reason="no direction was evaluated";
        std::vector<double> triedMasses;
        std::string massMove="accepted controller proposal";
        for(int proposal=0;proposal<3&&remaining()>0;++proposal)
        {
            triedMasses.push_back(currentMass);
            cb.SetCorrectionMass(currentMass);
            const bool cold=factorBase<0||newton-factorBase>o.maximumInverseAge;
            LinearResult result;
            Vec original;
            SteadyMeritRecord assessment;
            double priorForcing=1.;
            const std::array<double,4> forcingLevels={out.meritEnd->nextForcing,.1,.03,.01};
            for(double target:forcingLevels)
            {
                target=std::clamp(target,o.forcingMinimum,o.forcingMaximum);
                if(target>=priorForcing||remaining()<=0)continue;
                if(priorForcing<1.)
                {
                    if(std::max(assessment.linearNorm/norm,
                                assessment.scaledError/assessment.scaledBaseNorm)<=target)continue;
                    if((assessment.decision=="ACTION_UNRESOLVED"&&!assessment.modelEvaluated)||
                       assessment.actionUncertainty>=.5*target*norm||
                       assessment.scaledUncertainty>=.5*target*assessment.scaledBaseNorm)break;
                }
                const bool reuse=priorForcing<1.||!cold;
                const auto before=out.statistics.linearIterations;
                const Vec *warm=result.correction.size()?&result.correction:nullptr;
                result=Correct(x,f,norm,target,factorBase,reuse?"merit current inverse":
                    factorBase<0?"merit cold entry":"merit age cap",reuse,
                    std::min(remaining(),o.reuseLinearIterations),warm);
                assessment=MeasureMerit(x,f,norm,result,target,proposal,int(out.statistics.linearIterations-before),original);
                assessment.massMove=massMove;
                out.meritDirections.push_back(assessment);
                if(cb.ObserveMeritDirection)cb.ObserveMeritDirection(assessment);
                if(!result.converged&&priorForcing==1.&&reuse&&!refreshed&&remaining()>0&&
                   (factorBase!=newton||factorMass!=currentMass))
                {
                    refreshed=true;
                    const auto workBefore=out.statistics.linearIterations;
                    result=Correct(x,f,norm,target,factorBase,"bounded merit linear refresh",false,
                        std::min(remaining(),o.reuseLinearIterations),&result.correction);
                    assessment=MeasureMerit(x,f,norm,result,target,proposal,
                        int(out.statistics.linearIterations-workBefore),original);
                    assessment.massMove=massMove;
                    out.meritDirections.push_back(assessment);
                    if(cb.ObserveMeritDirection)cb.ObserveMeritDirection(assessment);
                }
                priorForcing=target;
                if(assessment.decision=="USEFUL_DESCENT")break;
            }
            reason=assessment.decision.empty()?"declared current linear work exhausted":assessment.decision;
            const double anchor=o.merit->reanchorMass;
            auto untried=[&](double mass) {
                return std::find(triedMasses.begin(),triedMasses.end(),mass)==triedMasses.end();
            };
            const bool declaredAnchor=anchor>currentMass&&untried(anchor);
            double alternative=declaredAnchor?anchor:StrongerMass(currentMass);
            if(!untried(alternative))alternative=.5*currentMass;
            const bool canChange=proposal<2&&remaining()>0&&untried(alternative);
            RecordLinear(result,priorForcing);
            if(reason=="USEFUL_DESCENT"||(reason=="WEAK_VERIFIED_DESCENT"&&!canChange))
            {
                const auto before=out.lineSearch.size();
                if(Backtrack(x,f,norm,PhysicalCorrection(result),original,assessment.step,
                             assessment.activeSlope,std::min(8,o.maximumLineSearchTrials),proposal))
                    return;
                const bool evaluated=std::any_of(out.lineSearch.begin()+before,out.lineSearch.end(),
                    [](const auto &trial){return trial.evaluated;});
                reason=evaluated?"FEASIBLE_ORIGINAL_LINE_SEARCH_EXHAUSTED":"NONLINEAR_DOMAIN_RESTRICTION";
            }
            if(!canChange)break;
            massMove=declaredAnchor?"declared re-anchor":
                alternative>currentMass?"bounded strengthening":"bounded weakening";
            currentMass=alternative;
        }
        throw Stop("ORIGINAL_MERIT_DIRECTION_STOP: "+reason);
    }

    bool Backtrack(Vec &x,Vec &f,double &norm,const Vec &correction,const Vec &originalAction,
                   double step,double modelSlope,int trialLimit,int proposal=-1)
    {
        if((o.measureDirection||o.merit)&&cb.InspectDirection)cb.InspectDirection(newton,x,correction);
        bool accepted=false;
        {
            Timer timer(out.statistics.lineSearchSeconds);
            for(int trial=0;trial<trialLimit;++trial)
            {
                Deadline();
                SteadyLineSearchRecord ls;
                ls.newton=newton;ls.trial=trial;ls.step=step;ls.proposal=proposal;
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
                    if(o.measureDirection||o.merit)
                    {
                        const Vec predicted=f+step*originalAction;
                        const Vec remainder=evaluation.residual-predicted;
                        ls.predictedNorm=Norm(predicted);
                        const double decrease=norm-*ls.predictedNorm;
                        if(decrease>0.)
                        {
                            const double ratio=(norm-ls.trueNorm)/decrease;
                            if(std::isfinite(ratio))ls.decreaseRatio=ratio;
                        }
                        if(cb.ResidualBlocks)
                        {
                            ls.predictedBlocks=cb.ResidualBlocks(predicted);
                            ls.actualBlocks=cb.ResidualBlocks(evaluation.residual);
                            ls.remainderBlocks=cb.ResidualBlocks(remainder);
                        }
                        if(cb.InspectTrial)cb.InspectTrial(newton,trial,step,candidate,remainder);
                    }
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
                    if(o.merit)AdvanceMerit(norm,ls);
                    x=candidate;f=std::move(evaluation.residual);norm=ls.trueNorm;
                    Observe(x,f,norm,true);accepted=true;
                    ++out.statistics.acceptedNewtonSteps;
                    auto &nr=out.newton.back();nr.accepted=true;
                    nr.acceptedNorm=norm;nr.acceptedStep=step;
                    nr.acceptedSeconds=Seconds(start);
                    break;
                }
                double reduction=o.lineSearchReduction;
                if(o.adaptiveTrialStep&&ls.evaluated&&modelSlope<0.)
                {
                    // Safeguarded minimizer of the quadratic through the
                    // measured trial with the current model slope at zero.
                    const double curvature=ls.trueNorm-norm-step*modelSlope;
                    if(curvature>0.)
                        reduction=std::clamp(-step*modelSlope/(2.*curvature),.1,.75);
                }
                step*=reduction;
            }
        }
        return accepted;
    }

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
        if(accepted)
        {
            if(out.lastAcceptedX.size()==0||!(x.array()==out.lastAcceptedX.array()).all())
                out.lastAcceptedNewtonIteration=newton;
            out.lastAcceptedX=x;
            out.lastAcceptedNorm=norm;
            if(cb.ObserveAccepted)cb.ObserveAccepted(x,norm);
        }
    }

    void Converged()
    {out.converged=true;out.status="CONVERGED";out.reason="current true steady residual meets caller acceptance norm";}

    SteadyKrylovResult Finish()
    {
        if(!out.newton.empty()&&!out.newton.back().linearEvaluated)
            out.newton.back().missingLinearReason+=": "+out.reason;
        out.statistics.totalSeconds=Seconds(start);return std::move(out);
    }

    LinearResult Correct(const Vec &x,const Vec &f,double norm,double forcing,
                         int &factorBase,const char *reason,bool reuse,int workLimit=-1,const Vec *guess=nullptr)
    {
        auto &record=out.inverses.emplace_back();
        record.base=newton;record.factorBase=reuse?factorBase:newton;
        record.age=newton-record.factorBase;record.reason=reason;
        record.rebuild=!reuse&&bool(cb.PreparePreconditioner);
        const auto before=out.statistics.linearIterations;
        try
        {
            if(record.rebuild)
            {
                out.newton.back().missingLinearReason="inverse preparation";
                ++out.statistics.preconditionerPreparations;
                Timer timer(out.statistics.preparationSeconds);
                cb.PreparePreconditioner(x,f);
                factorBase=newton;
                record.prepared=true;
                factorSerial=out.statistics.preconditionerPreparations;
                factorMass=currentMass;
            }
            out.newton.back().missingLinearReason="current linear solve";
            record.serial=factorSerial;
            if(o.merit)record.mass=factorMass;
            auto result=Linear(x,f,norm,forcing,workLimit>=0?workLimit:(reuse?o.reuseLinearIterations:0),guess);
            record.iterations=int(out.statistics.linearIterations-before);
            record.converged=result.converged;
            record.trueNorm=result.trueNorm;
            return result;
        }
        catch(...)
        {
            record.iterations=int(out.statistics.linearIterations-before);
            throw;
        }
    }

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

    LinearResult Linear(const Vec &x,const Vec &f,double baseNorm,double forcing,int workLimit,const Vec *guess)
    {
        Timer timer(out.statistics.linearSeconds);
        const Vec b=-(sr.array()*f.array()).matrix();
        CheckVector(b,sx.size(),"scaled Newton RHS");
        const double bNorm=b.stableNorm();
        if(!Positive(bNorm))throw Stop("LINEAR_SCALING: nonzero true residual has zero/nonfinite scaled norm");
        const int n=static_cast<int>(sx.size()),m=std::min(n,o.restartDimension);
        Vec best=guess?*guess:Vec::Zero(n).eval(),bestAction=Vec::Zero(n);
        double bestNorm=baseNorm,bestScaled=bNorm;
        int work=0;
        for(int cycle=0;cycle<o.maximumLinearCycles;++cycle)
        {
            Deadline();++out.statistics.linearCycles;
            const Vec baseCorrection=best;
            Vec baseAction=Vec::Zero(n);
            Vec r=b;
            if(baseCorrection.cwiseAbs().maxCoeff()!=0)
            {
                baseAction=Jacobian(x,f,baseCorrection,1.,"restart-residual",cycle,-1);
                r-=baseAction.cwiseProduct(sr);
            }
            CheckVector(r,n,"current scaled linear residual");
            const double beta=r.stableNorm();
            const Vec raw=(r.array()/sr.array()).matrix();
            const double initialTrue=Norm(raw);
            if(cycle==0&&guess)
            {
                bestNorm=initialTrue;
                bestScaled=beta;
                bestAction=baseAction;
            }
            SteadyLinearRecord restart;restart.newton=newton;restart.cycle=cycle;
            restart.iteration=-1;restart.trueChecked=true;restart.projectedNorm=beta;
            restart.trueScaledNorm=beta;restart.trueNorm=initialTrue;
            if(cb.ResidualBlocks)restart.trueBlocks=cb.ResidualBlocks(raw);
            restart.converged=beta<=forcing*bNorm&&initialTrue<=forcing*baseNorm;
            out.linear.push_back(restart);
            if(restart.converged)
                return {baseCorrection,true,initialTrue,std::move(baseAction)};
            if(!Positive(beta))throw Stop("LINEAR_SCALING: inconsistent zero/nonfinite Krylov residual");
            Mat v=Mat::Zero(n,m+1),z=Mat::Zero(n,m),h=Mat::Zero(m+1,m);
            v.col(0)=r/beta;
            const int available=workLimit?std::min(m,workLimit-work):m;
            for(int j=0;j<available;++j)
            {
                Deadline();++out.statistics.linearIterations;++work;
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
                if(breakdown||k==available||k%o.trueLinearCheckInterval==0||record.projectedNorm<=forcing*bNorm)
                {
                    const Vec candidate=baseCorrection+z.leftCols(k)*coefficients;
                    // This is a new nonlinear-residual FD evaluation in the
                    // accumulated direction, not the Arnoldi recurrence norm.
                    Vec action=Jacobian(x,f,candidate,1.,"true-linear-residual",cycle,j);
                    const Vec check=b-action.cwiseProduct(sr);
                    CheckVector(check,n,"recomputed linear residual");
                    record.trueChecked=true;record.trueScaledNorm=check.stableNorm();
                    const Vec rawCheck=(check.array()/sr.array()).matrix();
                    record.trueNorm=Norm(rawCheck);
                    if(cb.ResidualBlocks)record.trueBlocks=cb.ResidualBlocks(rawCheck);
                    record.converged=record.trueScaledNorm<=forcing*bNorm&&record.trueNorm<=forcing*baseNorm;
                    if(record.trueNorm<bestNorm||
                       (record.trueNorm==bestNorm&&record.trueScaledNorm<bestScaled))
                    {best=candidate;bestAction=action;bestNorm=record.trueNorm;bestScaled=record.trueScaledNorm;}
                    out.linear.push_back(record);Deadline();
                    if(record.converged)return {candidate,true,record.trueNorm,std::move(action)};
                }
                else out.linear.push_back(record);
                if(breakdown)break;
            }
            if(workLimit&&work>=workLimit)break;
        }
        return {best,false,bestNorm,std::move(bestAction)};
    }

    const SteadyKrylovCallbacks &cb;
    const Vec &sx,&sr;
    const SteadyKrylovOptions &o;
    Clock::time_point start;
    int newton=-1;
    long long factorSerial=0;
    double currentMass=0.,factorMass=0.;
    SteadyKrylovResult out;
};
}

std::string SteadyMeritState::ToJson()const
{
    namespace J=RANSWallSolveJson;
    return J::Object({{"schema",J::Quote(Identity)},{"used_mass",J::Number(usedMass)},
        {"next_mass",J::Number(nextMass)},{"reference_norm",J::Number(referenceNorm)},
        {"previous_norm",J::Number(previousNorm)},{"accepted_norm",J::Number(acceptedNorm)},
        {"used_forcing",J::Number(usedForcing)},{"next_forcing",J::Number(nextForcing)},
        {"step",J::Number(step)},{"actual_to_predicted_decrease",ratio?J::Number(*ratio):"null"},
        {"accepted",J::Number(accepted)},{"good_steps",J::Number(goodSteps)},
        {"decision",J::Quote(decision)}});
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
            o.trueLinearCheckInterval>0&&o.maximumLineSearchTrials>0&&o.maximumFDReductions>=0&&
            o.maximumInverseAge>=0&&o.reuseLinearIterations>0&&
            (o.maximumInverseAge==0||bool(m_callbacks.PreparePreconditioner))&&
            (!o.adaptiveTrialStep||o.measureDirection),
            "invalid declared steady algorithm cycle bounds");
    Require(Positive(o.fdRelativeStep)&&Positive(o.fdMinimumRelativeStep)&&
            o.fdMinimumRelativeStep<=o.fdRelativeStep&&Fraction(o.fdReduction)&&
            Fraction(o.lineSearchReduction)&&Fraction(o.armijo)&&o.armijo<.5,
            "invalid declared FD or norm-Armijo policy");
    Require(Fraction(o.forcingInitial)&&Fraction(o.forcingMinimum)&&Fraction(o.forcingMaximum)&&
            o.forcingMinimum<=o.forcingInitial&&o.forcingInitial<=o.forcingMaximum&&
            Positive(o.forcingGamma)&&Positive(o.forcingExponent)&&Fraction(o.arnoldiBreakdownRelative),
            "invalid declared inexact-Newton forcing or Arnoldi policy");
    if(o.merit)
    {
        const auto &m=*o.merit;
        Require(bool(m_callbacks.CorrectionAction)&&bool(m_callbacks.SetCorrectionMass)&&
                bool(m_callbacks.ResidualBlocks)&&!o.allowIncompleteLinearStep&&
                std::isfinite(m.initialMass)&&m.initialMass>=0.&&std::isfinite(m.reanchorMass)&&m.reanchorMass>=0.&&
                m.maximumWork>0,"invalid original-merit correction policy");
        if(m.restored)
        {
            const auto &r=*m.restored;
            Require(std::isfinite(r.usedMass)&&r.usedMass>=0.&&std::isfinite(r.nextMass)&&r.nextMass>=0.&&
                    Positive(r.referenceNorm)&&Positive(r.previousNorm)&&Positive(r.acceptedNorm)&&
                    r.usedForcing>=o.forcingMinimum&&r.usedForcing<=o.forcingMaximum&&
                    r.nextForcing>=o.forcingMinimum&&r.nextForcing<=o.forcingMaximum&&
                    r.accepted>=0&&r.goodSteps>=0&&r.goodSteps<2&&std::isfinite(r.step)&&r.step>=0.&&r.step<=1.&&
                    (!r.ratio||std::isfinite(*r.ratio))&&!r.decision.empty(),"invalid restored original-merit state");
        }
    }

}

SteadyKrylovResult SteadyKrylov::Solve(const SteadyVector &initial)const
{return Runner(m_callbacks,m_unknownScale,m_residualRowScale,m_options).Solve(initial);}

SteadyDerivativeProbe SteadyKrylov::ProbeDirection(const SteadyVector &base,
    const SteadyVector &direction,const std::vector<double> &multipliers)const
{return Runner(m_callbacks,m_unknownScale,m_residualRowScale,m_options).Probe(base,direction,multipliers);}
}
