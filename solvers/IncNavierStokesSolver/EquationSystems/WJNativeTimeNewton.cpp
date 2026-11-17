// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeCoupling.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeFailures.h>
#include <chrono>
#include <cmath>
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
    catch(...){save();throw;}

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
            { AccountTime timer(info.timing.rebuild);ResetPreconditioners(); }
            CoupledState initial;
            { AccountTime timer(info.timing.wall);initial=DirichletState(history.current,time); }
            CoupledState pcState;
            Sparse epsilonPCFull;
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
                Sparse flowPC=VelocityTangent(data,info.timing), kPC, ePC;
                {
                    AccountTime timer(info.timing.assembly);
                    kPC=m_k.Diffusion(data.kRules,data.Dk);
                    ePC=m_e.Diffusion(m_rules.Static(),data.De,data.reactionE);
                    epsilonPCFull=ePC+alpha*m_e.Mass();
                }
                {
                    AccountTime timer(info.timing.rebuild);
                    m_flow.PrepareLinearPreconditioner(m_flow.System(alpha,flowPC,1.));++info.rebuilds;
                    m_k.PrepareLinearPreconditioner(m_k.System(alpha,kPC));++info.rebuilds;
                    m_e.PrepareLinearPreconditioner(m_e.System(alpha,ePC));++info.rebuilds;
                }
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
                    Vec cross=Vec::Zero(2*m_flow.VelocitySize()+m_flow.PressureSize());
                    cross.head(2*m_flow.VelocitySize())=m_flow.KnownScalarKLoad(dk);
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
                    result.tail(ne)=m_e.ApplyLinearPreconditioner(
                        rhs.tail(ne)-m_e.Restrict(epsilonPCFull*trace));
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
                nextAdvection=Advection(state,final.kRules);
            }
            auto next=history;
            next.previous=history.current;next.current=state;
            next.previousAdvection=history.advection;next.advection=std::move(nextAdvection);
            ++next.step;next.time=time;
            CheckHistory(next);
            checkDeadline();
            history=std::move(next);
            save();
            return info;
        }
        catch(const NativeLinearFailure &error)
        {
            info.retryReason=error.what();
            if(attempt==1){save();throw;}
        }
        catch(...){save();throw;}
    }
    throw std::runtime_error("unreachable implicit retry state");
}
}
