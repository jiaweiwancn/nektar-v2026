// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeCoupling.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeFailures.h>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <stdexcept>

namespace Nektar::WJNative
{
namespace
{
double Seconds()
{
    return std::chrono::duration<double>(
        std::chrono::steady_clock::now().time_since_epoch()).count();
}
struct AccountTime
{
    explicit AccountTime(double &destination) : value(destination), start(Seconds()) {}
    ~AccountTime() { value += Seconds()-start; }
    double &value;
    double start;
};
constexpr std::array<const char*,7> RefreshNames{
    "cold","alpha","drift","age","control","retry","injected"};
double OperatorChange(const Sparse &current,const Sparse &saved)
{
    if(current.rows()!=saved.rows()||current.cols()!=saved.cols())
        throw std::logic_error("Newton inverse-operator cache has different field maps");
    const double scale=saved.norm(), difference=(current-saved).norm();
    if(!(std::isfinite(scale)&&scale>0&&std::isfinite(difference)&&difference>=0))
        throw std::runtime_error("invalid Newton inverse-operator drift norm");
    const double relative=difference/scale;
    if(!std::isfinite(relative))throw std::runtime_error("nonfinite Newton inverse-operator drift");
    return relative;
}
double DataChange(const Vec &current,const Vec &saved)
{
    if(current.size()!=saved.size()||current.size()==0)
        throw std::logic_error("Newton inverse cache has different current-data inventories");
    const double scale=saved.norm(),difference=(current-saved).norm();
    if(!(std::isfinite(scale)&&scale>0&&std::isfinite(difference)&&difference>=0))
        throw std::runtime_error("invalid current-data inverse-cache drift norm");
    const double relative=difference/scale;
    if(!std::isfinite(relative))throw std::runtime_error("nonfinite current-data inverse-cache drift");
    return relative;
}
}

Vec Coupling::TimePreconditionerData(const Data &data)const
{
    const auto&p=m_reference->Parameters();
    const double K0=p.KIn,A0=m_reference->Amplitude(K0),c0=m_reference->Cutoff(K0);
    const auto&rules=m_rules.Static();
    if(data.kCoefficients.size()==0||data.velocityInputs.size()!=rules.size()||data.coefficients.size()!=2*data.kCoefficients.size())
        throw std::logic_error("current-data inverse-cache static point inventory differs");
    Vec value(10*data.kCoefficients.size()+data.eCoefficients.size());
    Eigen::Index point=0,epsilon=0,entry=0;
    for(std::size_t i=0;i<rules.size();++i)
    {
        const auto&rule=*rules[i];const bool outer=m_e.Owns(rule.geometry);
        if(data.velocityInputs[i].size()!=std::size_t(rule.w.size()))
            throw std::logic_error("current-data inverse-cache owned quadrature differs");
        for(const auto&input:data.velocityInputs[i])
        {
            if(point>=data.kCoefficients.size()||entry+10+2*int(outer)>value.size())
                throw std::logic_error("current-data inverse-cache tuple length differs");
            const double E0=outer?A0/input.r:A0/std::max(input.r,c0);
            const double Dk0=p.nu+p.Cd*K0*K0/(p.sigmaK*E0);
            const double De0=p.nu+p.Cd*K0*K0/(p.sigmaE*E0);
            const double reaction0=1.82*E0/K0;
            if(!(std::isfinite(E0)&&E0>0&&std::isfinite(Dk0)&&Dk0>0&&
                 std::isfinite(De0)&&De0>0&&std::isfinite(reaction0)&&reaction0>0))
                throw std::runtime_error("invalid fixed current-data inverse-cache units");
            value[entry++]=input.k/K0;value[entry++]=input.E/E0;value[entry++]=input.aux/Dk0;
            for(int component:{0,2})for(int derivative:{0,2})
                value[entry++]=input.gradient[component][derivative]*K0/E0;
            value[entry++]=data.coefficients[2*point]/Dk0;
            value[entry++]=data.coefficients[2*point+1]/Dk0;
            value[entry++]=data.kCoefficients[point++]/Dk0;
            if(outer)
            {
                if(epsilon+2>data.eCoefficients.size())
                    throw std::logic_error("current-data inverse-cache outer inventory differs");
                value[entry++]=data.eCoefficients[epsilon++]/De0;
                value[entry++]=data.eCoefficients[epsilon++]/reaction0;
            }
        }
    }
    if(entry!=value.size()||point!=data.kCoefficients.size()||epsilon!=data.eCoefficients.size()||!value.allFinite())
        throw std::logic_error("incomplete/nonfinite current-data inverse-cache tuple");
    if(!(std::isfinite(value.norm())&&value.norm()>0))
        throw std::runtime_error("nonpositive/nonfinite current-data inverse-cache tuple norm");
    return value;
}

void WriteStepPreconditionerFields(std::ostream &out,const CoupledOptions &options,
                                  const CoupledStepInfo &info)
{
    const bool amdOrdering=options.stepPreconditioner=="CurrentDataAMDOrderingV1";
    out<<std::setprecision(17)<<",\"preconditioner_policy\":\""<<options.stepPreconditioner
       <<"\",\"mixed_sparse_ordering\":\""<<(amdOrdering?"AMD":"COLAMD")
       <<"\",\"drift_measure\":\""<<(amdOrdering||options.stepPreconditioner=="CurrentDataDriftV1"?
            "ScaledCurrentQuadratureDataL2V1":"CurrentInverseOperatorFrobeniusV1")
       <<"\",\"relative_drift_limit\":"<<options.stepPCRelativeDrift
       <<",\"maximum_uses\":"<<options.stepPCMaximumUses
       <<",\"refreshes\":"<<info.timePCRefreshes<<",\"reuses\":"<<info.timePCReuses
       <<",\"maximum_used_age\":"<<info.timePCMaximumUsedAge<<",\"refresh_reasons\":{";
    for(int i=0;i<int(RefreshNames.size());++i)
    {if(i)out<<',';out<<'"'<<RefreshNames[i]<<"\":"<<info.timePCRefreshReasons[i];}
    out<<"},\"preconditioner_events_recorded\":"<<(options.recordNonlinearTrace&&options.advancement=="CoupledJFNKV1"?"true":"false")
       <<",\"preconditioner_events\":[";
    for(int i=0;i<int(info.timePCEvents.size());++i)
    {
        const auto &e=info.timePCEvents[i];if(i)out<<',';
        out<<"{\"attempt\":"<<e.attempt<<",\"base\":"<<e.base
           <<",\"alpha\":"<<e.alpha<<",\"drift\":"<<e.drift
           <<",\"age_before\":"<<e.ageBefore<<",\"age_used\":"<<e.ageUsed
           <<",\"reason\":\""<<e.reason<<"\",\"refreshed\":"<<(e.refreshed?"true":"false")<<'}';
    }
    out<<']';
}

CoupledStepInfo Coupling::AdvanceNewton(CoupledHistory &history, bool injectRetry)
{
    const double began=Seconds();
    CheckHistory(history);
    m_constitutiveStats={};
    m_lastInputHistory=history;
    m_hasInputHistory=true;
    m_lastAttempt={};
    CoupledStepInfo info;
    const double time=(history.step+1)*history.dt;
    const double alpha=(history.step==0?1.:1.5)/history.dt;
    const int nf=m_flow.Size(), nk=m_k.FreeSize(), ne=m_e.FreeSize();
    const int size=nf+nk+ne;
    auto save=[&] {
        info.seconds=Seconds()-began;
        info.constitutiveStats=m_constitutiveStats;
        info.timing.other=info.seconds-info.timing.Sum();
        m_lastAttempt=info;
    };
    const double deadline=began+m_options.stepNewton.maximumSeconds;
    auto checkDeadline=[&] {
        if(Seconds()>=deadline)
            throw std::runtime_error("implicit physical step exhausted cooperative deadline; accepted history preserved");
    };
    Loads external;
    Vec unknownScale, rowScale(size);
    try
    {
        external=Forcing(time,info.timing);
        AccountTime timer(info.timing.nonlinear);
        const double currentEX=history.step==0?1.:2.;
        const double previousEX=history.step==0?0.:-1.;
        external.flow-=currentEX*history.advection.flow+previousEX*history.previousAdvection.flow;
        external.k-=currentEX*history.advection.k+previousEX*history.previousAdvection.k;
        external.epsilon-=currentEX*history.advection.epsilon+previousEX*history.previousAdvection.epsilon;
        unknownScale=UnknownScales();
        rowScale.head(nf)=m_flow.ResidualRowScale();
        rowScale.segment(nf,nk)=m_k.ResidualRowScale();
        rowScale.tail(ne)=m_e.ResidualRowScale();
        checkDeadline();
    }
    catch(...){
        {AccountTime timer(info.timing.rebuild);ResetPreconditioners();}
        save();throw;
    }

    auto norms=[&](const Vec &raw) {
        if(raw.size()!=size||!raw.allFinite())
            throw std::runtime_error("invalid current implicit Newton residual");
        const auto flow=m_flow.RawResidualNorm(raw.head(nf));
        return CoupledResidual{flow.momentum,flow.continuity,flow.gauge,
            m_k.FreeResidualNorm(raw.segment(nf,nk)),
            m_e.FreeResidualNorm(raw.tail(ne)),0};
    };

    // Only a typed linear failure receives the same-dt cold retry. Neither
    // Newton iteration limits nor scientific/final guard failures are retries.
    for(int attempt=0;attempt<2;++attempt)
    {
        if(attempt)++info.retries;
        try
        {
            checkDeadline();
            if(attempt||m_options.stepPreconditioner=="RebuildEveryNewtonV1")
            { AccountTime timer(info.timing.rebuild);ResetPreconditioners(); }
            CoupledState initial;
            { AccountTime timer(info.timing.wall);initial=DirichletState(history.current,time); }
            CoupledState pcState;
            const Sparse *epsilonPCFull=nullptr;
            int prepareBase=0;
            Vec validatedUnknowns;
            std::shared_ptr<const Data> validatedData;
            auto unpack=[&](const Vec &x) {
                AccountTime timer(info.timing.wall);
                return FromFree(x,time);
            };
            auto currentData=[&](const Vec &x,const CoupledState &s) {
                if(validatedData&&x.size()==validatedUnknowns.size()&&
                   (x.array()==validatedUnknowns.array()).all())return validatedData;
                return std::shared_ptr<const Data>(
                    std::make_shared<Data>(Evaluate(s,time,false,info.timing)));
            };
            auto equations=[&](const Vec &x) {
                const double evaluationBegan=Seconds();
                const auto state=unpack(x);
                const auto owned=currentData(x,state);
                AccountTime timer(info.timing.residual);
                auto action=CurrentActions(state,*owned);
                // Stable coefficient differences include the complete
                // essential/Gamma history before the physical mass action.
                CoupledState delta{alpha*(state.flow-history.current.flow),
                    alpha*(state.k-history.current.k),
                    alpha*(state.epsilon-history.current.epsilon)};
                if(history.step>0)
                {
                    const double a=.5/history.dt;
                    delta.flow+=a*(history.previous.flow-history.current.flow);
                    delta.k+=a*(history.previous.k-history.current.k);
                    delta.epsilon+=a*(history.previous.epsilon-history.current.epsilon);
                }
                action.flow+=m_flow.Core().MassAction(delta.flow)-external.flow;
                action.k+=m_k.Mass()*delta.k-external.k;
                action.epsilon+=m_e.Mass()*delta.epsilon-external.epsilon;
                Vec raw(size);
                raw.head(nf)=m_flow.Restrict(action.flow);
                if(m_flow.Gauge())
                {
                    raw.segment(m_flow.PressureOffset(),m_flow.PressureSize())+=
                        m_flow.PressureMean()*x[nf-1];
                    raw[nf-1]=m_flow.PressureMean().dot(
                        state.flow.tail(m_flow.PressureSize()));
                }
                raw.segment(nf,nk)=m_k.Restrict(action.k);
                raw.tail(ne)=m_e.Restrict(action.epsilon);
                if(m_options.recordNonlinearTrace&&validatedData&&
                   (x.array()==validatedUnknowns.array()).all())
                    info.nonlinearTrace.push_back({norms(raw),false,Seconds()-evaluationBegan});
                return raw;
            };

            SteadyKrylovCallbacks callbacks;
            callbacks.ValidTrial=[&](const Vec &x) {
                try
                {
                    const auto state=unpack(x);
                    AccountTime timer(info.timing.guards);
                    const auto k=m_k.PolynomialBounds(state.k);
                    const auto e=m_e.PolynomialBounds(state.epsilon);
                    if(k[0]<=0||e[0]<=0)
                        return SteadyTrialValidity{false,"nonpositive transported scalar polynomial"};
                    return SteadyTrialValidity{};
                }
                catch(const NativeInadmissibleState &error)
                { return SteadyTrialValidity{false,error.what()}; }
            };
            // These states are uncommitted iterations, exactly as in Picard.
            // Root/domain/covariance are current; a physical acceptance below
            // always recomputes the complete principal certificate afresh.
            callbacks.ValidateNonlinearState=[&](const Vec &x) {
                try
                {
                    auto full=std::make_shared<Data>(Evaluate(unpack(x),time,false,info.timing));
                    validatedUnknowns=x;
                    validatedData=std::move(full);
                }
                catch(const NativeInadmissibleState &error)
                { throw SteadyInadmissibleTrial(error.what()); }
            };
            callbacks.Residual=[&](const Vec &x)->Vec {
                ++info.coupledResidualEvaluations;
                try { return equations(x); }
                catch(const NativeInadmissibleState &error)
                { throw SteadyInadmissibleTrial(error.what()); }
            };
            callbacks.TrueNorm=[&](const Vec &raw) {
                AccountTime timer(info.timing.residual);
                return norms(raw).Maximum();
            };
            callbacks.PreparePreconditioner=[&](const Vec &x,const Vec &) {
                pcState=unpack(x);
                const auto owned=currentData(x,pcState);
                const auto &data=*owned;
                // EX advection is fixed history, so it has no current-state
                // transport Jacobian. Alpha is the actual BE/BDF2 coefficient.
                Sparse flowPC,kPC,ePC;Vec currentTuple;
                const bool amdOrdering=m_options.stepPreconditioner=="CurrentDataAMDOrderingV1";
                const bool currentDataPolicy=amdOrdering||m_options.stepPreconditioner=="CurrentDataDriftV1";
                auto assembleInverseOperators=[&] {
                    flowPC=VelocityTangent(data,info.timing);
                    AccountTime timer(info.timing.assembly);
                    kPC=m_k.Diffusion(data.kRules,data.Dk);
                    ePC=m_e.Diffusion(m_rules.Static(),data.De,data.reactionE);
                };
                double drift=0;
                const int previousAge=m_timePCUses;
                if(currentDataPolicy)
                {
                    AccountTime timer(info.timing.assembly);
                    currentTuple=TimePreconditionerData(data);
                    if(m_timePCReady)drift=DataChange(currentTuple,m_timePCData);
                }
                else
                {
                    assembleInverseOperators();
                    AccountTime timer(info.timing.assembly);
                    if(m_timePCReady)
                        drift=std::max({OperatorChange(flowPC,m_timePCOperators[0]),
                            OperatorChange(kPC,m_timePCOperators[1]),
                            OperatorChange(ePC,m_timePCOperators[2])});
                }
                info.coefficientDrift=std::max(info.coefficientDrift,drift);
                // This norm is a numerical refresh heuristic for inverse
                // operators only. No residual, Jv, load or physical coefficient
                // is ever taken from this cache. Alpha changes always refresh.
                int reason=-1;
                if(injectRetry&&attempt==0)reason=6;
                else if(attempt)reason=5; // cold retry: fresh at EVERY base
                else if(m_options.stepPreconditioner=="RebuildEveryNewtonV1")reason=4;
                else if(!m_timePCReady)reason=0;
                else if(alpha!=m_timePCAlpha)reason=1;
                else if(drift>m_options.stepPCRelativeDrift)reason=2;
                else if(m_timePCUses>=m_options.stepPCMaximumUses)reason=3;
                if(reason>=0)
                {
                    m_timePCReady=false; // a partial factorization cannot be reused
                    if(currentDataPolicy)assembleInverseOperators();
                    {
                        AccountTime timer(info.timing.rebuild);
                        m_flow.PrepareLinearPreconditioner(m_flow.System(alpha,flowPC,1.),amdOrdering);++info.rebuilds;
                        m_k.PrepareLinearPreconditioner(m_k.System(alpha,kPC));++info.rebuilds;
                        m_e.PrepareLinearPreconditioner(m_e.System(alpha,ePC));++info.rebuilds;
                    }
                    {
                        AccountTime timer(info.timing.assembly);
                        m_timePCOperators={std::move(flowPC),std::move(kPC),std::move(ePC)};
                        m_timePCEpsilonFull=m_timePCOperators[2]+alpha*m_e.Mass();
                        m_timePCData=std::move(currentTuple);
                    }
                    m_timePCAlpha=alpha;m_timePCUses=0;m_timePCReady=true;
                    ++info.timePCRefreshes;++info.timePCRefreshReasons[reason];
                }
                else ++info.timePCReuses;
                ++m_timePCUses;
                info.timePCMaximumUsedAge=std::max(info.timePCMaximumUsedAge,m_timePCUses);
                epsilonPCFull=&m_timePCEpsilonFull;
                if(m_options.recordNonlinearTrace)
                    info.timePCEvents.push_back({attempt,prepareBase,previousAge,m_timePCUses,
                        alpha,drift,reason<0?"reuse":RefreshNames[reason],reason>=0});
                ++prepareBase;
                if(injectRetry&&attempt==0)
                    throw NativeLinearFailure("qualification injected recoverable Krylov failure");
            };
            callbacks.Precondition=[&](const Vec &rhs) {
                ++info.coupledPreconditionerApplications;
                Vec result(size), dk;
                {
                    AccountTime timer(info.timing.scalar);
                    result.segment(nf,nk)=m_k.ApplyLinearPreconditioner(rhs.segment(nf,nk));
                    dk=m_k.Expand(result.segment(nf,nk),Vec::Zero(m_k.Size()));
                }
                {
                    AccountTime timer(info.timing.flow);
                    Vec cross;
                    {
                        AccountTime coupling(info.timingDetail.scalarMomentumCoupling);
                        cross = Vec::Zero(2*m_flow.VelocitySize()+m_flow.PressureSize());
                        cross.head(2*m_flow.VelocitySize()) = m_flow.KnownScalarKLoad(dk);
                    }
                    result.head(nf)=m_flow.ApplyLinearPreconditioner(rhs.head(nf)+m_flow.Restrict(cross));
                }
                Vec trace;
                {
                    AccountTime timer(info.timing.wall);
                    trace=m_e.EssentialLift([&](double x,double,const Vector &normal) {
                        if(normal[2]>=-.5)return 0.;
                        const Element *owner=nullptr;
                        for(const auto &e:m_flow.Elements())
                            if(std::abs(e.rl-m_geometry.H)<1e-12&&x>=e.xl-1e-12&&x<=e.xr+1e-12)
                            { owner=&e;break; }
                        if(!owner)throw std::runtime_error("implicit Gamma derivative has no outer owner");
                        const auto rule=m_rules.Points(*owner,Vec::Constant(1,x),
                            Vec::Constant(1,m_geometry.d0+m_geometry.H));
                        return m_reference->AmplitudeDerivative(InterfaceK(*owner,x,pcState))*
                            m_k.Evaluate(rule,dk)[0]/m_geometry.H;
                    },m_rules);
                }
                {
                    AccountTime timer(info.timing.scalar);
                    if(!epsilonPCFull||!m_timePCReady)
                        throw std::logic_error("unprepared Newton inverse-operator cache");
                    result.tail(ne)=m_e.ApplyLinearPreconditioner(
                        rhs.tail(ne)-m_e.Restrict((*epsilonPCFull)*trace));
                }
                return result;
            };
            checkDeadline();
            auto policy=m_options.stepNewton;
            policy.maximumSeconds=deadline-Seconds();
            SteadyKrylov solver(std::move(callbacks),unknownScale,rowScale,policy);
            SteadyKrylovResult result;
            const double engineBegan=Seconds(), accountedBefore=info.timing.Sum();
            auto accountEngine=[&] {
                // Exclusive Newton/GMRES/Armijo orchestration: callbacks
                // already own their disjoint equation/closure/guard work.
                info.timing.nonlinear+=std::max(0.,Seconds()-engineBegan-
                    (info.timing.Sum()-accountedBefore));
            };
            try { result=solver.Solve(PackFree(initial)); }
            catch(...){accountEngine();info.coupledStatisticsComplete=false;throw;}
            accountEngine();
            const auto &stats=result.statistics;
            info.nonlinearIterations+=int(stats.acceptedNewtonSteps);
            info.coupledKrylovIterations+=stats.linearIterations;
            info.coupledJacobianActions+=stats.jacobianActions;
            info.coupledRejectedTrials+=stats.fdRejectedTrials+stats.lineSearchRejectedTrials;
            if(result.hasEvaluatedState)info.residual=norms(result.residual);
            if(!result.converged)
            {
                if(result.reason.rfind("LINEAR_NOT_CONVERGED:",0)==0)
                    throw NativeLinearFailure(result.reason);
                throw std::runtime_error("current implicit Newton step failed: "+result.reason);
            }
            const auto state=unpack(result.x);
            checkDeadline();
            const auto final=Evaluate(state,time,true,info.timing);
            info.residual=Residual(state,final,external,alpha,time,info.timing,&history);
            if(info.residual.Maximum()>m_options.nonlinearTolerance)
                throw std::runtime_error("fresh current implicit residual failed after eliminating gauge multiplier");
            checkDeadline();
            info.minimumK=final.minK;info.minimumEpsilon=final.minE;
            info.minimumCovarianceOverK=final.minCovariance;info.minimumPrincipalOverMu=final.minPrincipal;
            CoupledState nextAdvection;
            {
                AccountTime timer(info.timing.assembly);
                nextAdvection = Advection(state, final.kRules, &info.timingDetail);
            }
            {
                AccountTime movement(info.timingDetail.historyMovement);
                auto next=history;
                next.previous=history.current;next.current=state;
                next.previousAdvection=history.advection;next.advection=std::move(nextAdvection);
                ++next.step;next.time=time;
                CheckHistory(next);
                checkDeadline();
                history=std::move(next);
            }
            save();
            return info;
        }
        catch(const NativeLinearFailure &error)
        {
            info.retryReason=error.what();
            {AccountTime timer(info.timing.rebuild);ResetPreconditioners();}
            if(attempt==1){save();throw;}
        }
        catch(...){
            {AccountTime timer(info.timing.rebuild);ResetPreconditioners();}
            save();throw;
        }
    }
    throw std::runtime_error("unreachable implicit retry state");
}
}
