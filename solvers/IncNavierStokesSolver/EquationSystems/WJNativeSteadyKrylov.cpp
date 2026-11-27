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

struct AuxiliaryRhs
{
    const Vec &defect;
    const std::array<double,6> &blockBudgets;
    double scaledBudget;
    int maximumProducts;
};

struct CurvedCandidate
{
    Vec x;
    Evaluation evaluation;
};

using FDSample=std::pair<double,Vec>;

bool FDStop(const Stop &error)
{return std::string(error.what()).rfind("FD_",0)==0;}

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
            if(o.frozenSizeColumns)
            {
                newton=0;
                SteadyNewtonRecord record;
                record.baseNorm=norm;
                record.forcing=o.forcingInitial;
                out.newton.push_back(record);
                if(cb.PrepareCorrection)cb.PrepareCorrection(x,f);
                int factorBase=-1;
                Correct(x,f,norm,o.forcingInitial,factorBase,"frozen size basis",false,o.frozenSizeColumns);
                out.status="FROZEN_SIZE_BASIS";
                out.reason="one frozen correction basis collected";
                return Finish();
            }
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
                if(o.merit&&o.merit->bounded)
                    BoundedDirections(x,f,norm,factorBase);
                else if(o.merit)
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

    void Audit(const Vec &x,const Vec &direction)
    {
        CheckVector(x,sx.size(),"audit base");
        CheckVector(direction,sx.size(),"physical audit direction");
        Require(direction.cwiseAbs().maxCoeff()>0.,"audit requires a nonzero direction");
        const Vec f=Base(x);
        if(cb.PrepareCorrection)cb.PrepareCorrection(x,f);
        AuditSpacings("seed-direction",x,f,direction.cwiseQuotient(sx),-1,-1,0);
    }

    void AuditDirections(const Vec &x,const Mat &directions)
    {
        CheckVector(x,sx.size(),"frozen candidate base");
        const Vec f=Base(x);
        if(cb.PrepareCorrection)cb.PrepareCorrection(x,f);
        newton=0;
        const double physicalTarget=o.forcingInitial*Norm(f);
        const double scaledTarget=o.forcingInitial*sr.cwiseProduct(f).stableNorm();
        for(int index=0;index<directions.cols();++index)
        {
            const Vec v=directions.col(index).cwiseQuotient(sx);
            auto audit=AuditData("frozen-size-candidate",x,f,v,-1,index);
            audit.proposal=index;
            audit.measurements={{"physical_target",physicalTarget},{"scaled_target",scaledTarget},
                {"scaled_size",v.cwiseAbs().maxCoeff()},{"forcing",o.forcingInitial}};
            const auto first=out.finiteDifferences.size();
            Vec previous, primary, half;
            double primaryStep=0.,halfStep=0.;
            const std::array<double,3> multipliers{1.,.5,2.};
            for(int sample=0;sample<3;++sample)
            {
                const auto label=std::to_string(sample);
                try
                {
                    const Vec action=Jacobian(x,f,v,multipliers[sample],"candidate-"+label,-1,index);
                    if(o.frozenFiniteTrialLimit&&action.size()&&out.finiteDifferences.size()>first)
                    {
                        if(sample==0){primary=action;primaryStep=out.finiteDifferences.back().parameterStep;}
                        if(sample==1){half=action;halfStep=out.finiteDifferences.back().parameterStep;}
                    }
                    audit.arrays.emplace_back("raw_action_"+label,action);
                    audit.arrays.emplace_back("raw_model_residual_"+label,f+action);
                    if(previous.size())
                        audit.arrays.emplace_back("raw_spacing_change_"+label,action-previous);
                    previous=action;
                }
                catch(const Stop &error)
                {
                    if(!FDStop(error))throw;
                    audit.notes.emplace_back("unavailable_sample_"+label,error.what());
                }
            }
            audit.stencils.assign(out.finiteDifferences.begin()+first,out.finiteDifferences.end());
            cb.ObserveActionAudit(audit);
            if(o.frozenFiniteTrialLimit)
                AuditFiniteTrials(x,f,directions.col(index),primary,half,primaryStep,halfStep,index);
        }
    }

    void AuditColumn(const Vec &x,const Vec &direction,const SteadyColoredColumnProbe &column,
                     const std::string &label)
    {
        const Vec f=Base(x);
        if(cb.PrepareCorrection)cb.PrepareCorrection(x,f);
        const auto started=Clock::now();
        const auto displacement=column.colorDirection.cwiseQuotient(sx).eval();
        Vec unit=Vec::Zero(sx.size());
        unit[column.column]=1.;
        const Vec single=Jacobian(x,f,unit,1.,"singleton",-1,-1);
        const Vec colored=Jacobian(x,f,displacement,1.,"complete-color",-1,-1);
        Require(out.finiteDifferences.size()==2&&
                out.finiteDifferences[0].method=="central"&&out.finiteDifferences[1].method=="central",
                "selected color comparison requires two central stencils at the same prescribed spacing");
        Vec outside=single-column.exactSingleAction, owned=Vec::Zero(sx.size());
        for(int row:column.support)
        {
            owned[row]=colored[row]-column.exactColorAction[row]-outside[row];
            outside[row]=0.;
        }
        auto data=AuditData(label.c_str(),x,f,unit,-1,-1);
        const double weight=direction[column.column]/sx[column.column];
        data.measurements={{"column",double(column.column)},{"color_columns",double(column.colorColumns)},
                           {"owned_rows",double(column.support.size())},{"coefficient_weight",weight},
                           {"current_fd",o.fdRelativeStep},{"audit_seconds",Seconds(started)}};
        data.arrays.emplace_back("raw_full_row_singleton",single);
        data.arrays.emplace_back("raw_full_row_color",colored);
        data.arrays.emplace_back("color_direction",column.colorDirection);
        data.arrays.emplace_back("raw_exact_singleton",column.exactSingleAction);
        data.arrays.emplace_back("raw_exact_color",column.exactColorAction);
        data.arrays.emplace_back("raw_outside_support",outside);
        data.arrays.emplace_back("raw_owned_color_delta",owned);
        data.arrays.emplace_back("raw_weighted_outside_support",weight*outside);
        data.arrays.emplace_back("raw_weighted_owned_color_delta",weight*owned);
        data.stencils=out.finiteDifferences;
        cb.ObserveActionAudit(data);
    }

private:

    SteadyKrylovActionAudit AuditData(const char *stage,const Vec &x,const Vec &f,
                                const Vec &direction,int cycle,int iteration) const
    {
        SteadyKrylovActionAudit data;
        data.stage=stage;data.newton=newton;data.cycle=cycle;data.iteration=iteration;
        data.factorSerial=factorSerial;
        data.arrays={{"base",x},{"raw_base_residual",f},{"unknown_scale",sx},
                     {"row_scale",sr},{"direction",sx.cwiseProduct(direction)}};
        return data;
    }

    void AuditFiniteTrials(const Vec &x,const Vec &f,const Vec &direction,const Vec &primary,
                           const Vec &half,double primaryStep,double halfStep,int proposal)
    {
        const double base=Norm(f), scaledBase=sr.cwiseProduct(f).stableNorm();
        const double physicalTarget=o.forcingInitial*base,scaledTarget=o.forcingInitial*scaledBase;
        const bool haveAllowance=primaryStep>0.&&halfStep>0.&&halfStep<.8*primaryStep;
        Vec allowance;
        double physicalAllowance=0.,scaledAllowance=0.;
        if(haveAllowance)
        {
            allowance=2.*(half-primary);
            physicalAllowance=Norm(allowance);
            scaledAllowance=sr.cwiseProduct(allowance).stableNorm();
        }
        const bool resolved=haveAllowance&&physicalAllowance<.5*physicalTarget&&scaledAllowance<.5*scaledTarget;
        for(int i=0;i<o.frozenFiniteTrialLimit;++i)
        {
            Deadline();
            const double fraction=std::ldexp(1.,-i);
            const Vec candidate=x+fraction*direction;
            const Vec displacement=(candidate-x).cwiseQuotient(sx);
            auto audit=AuditData("frozen-size-finite-trial",x,f,direction.cwiseQuotient(sx),-1,i);
            audit.proposal=proposal;
            audit.measurements={{"fraction",fraction},{"scaled_size",displacement.cwiseAbs().maxCoeff()},
                {"physical_target",physicalTarget},{"scaled_target",scaledTarget},{"armijo_limit",(1.-o.armijo*fraction)*base},
                {"physical_base",base},{"scaled_base",scaledBase},{"action_resolved",resolved?1.:0.}};
            audit.arrays.emplace_back("trial_free",candidate);
            audit.arrays.emplace_back("canonical_input_direction",direction);
            if(cb.ScaledUnknownBlocks)
            {
                const auto blocks=cb.ScaledUnknownBlocks(displacement);
                audit.arrays.emplace_back("scaled_displacement_blocks",Eigen::Map<const Vec>(blocks.data(),blocks.size()));
            }
            Vec prediction;
            double predicted=0.,scaledPredicted=0.;
            if(primary.size())
            {
                prediction=f+fraction*primary;
                predicted=Norm(prediction);
                scaledPredicted=sr.cwiseProduct(prediction).stableNorm();
                audit.arrays.emplace_back("raw_model_residual",prediction);
                audit.measurements.emplace_back("physical_model",predicted);
                audit.measurements.emplace_back("scaled_model",scaledPredicted);
            }
            if(haveAllowance)
            {
                audit.arrays.emplace_back("raw_action_allowance",fraction*allowance);
                audit.measurements.emplace_back("physical_action_allowance",fraction*physicalAllowance);
                audit.measurements.emplace_back("scaled_action_allowance",fraction*scaledAllowance);
            }
            auto evaluation=Trial(candidate,false);
            audit.measurements.emplace_back("nonlinear_validated",evaluation.nonlinearValidated?1.:0.);
            std::string status="FINITE_STATE_REJECTED";
            if(!evaluation.admissible)
                audit.notes.emplace_back("failure",evaluation.reason);
            bool useful=false;
            if(evaluation.admissible)
            {
                const double actual=Norm(evaluation.residual),scaledActual=sr.cwiseProduct(evaluation.residual).stableNorm();
                const bool armijo=actual<=(1.-o.armijo*fraction)*base;
                const bool root=actual<=o.acceptanceNorm;
                audit.arrays.emplace_back("raw_actual_residual",evaluation.residual);
                audit.measurements.emplace_back("physical_actual",actual);
                audit.measurements.emplace_back("scaled_actual",scaledActual);
                audit.measurements.emplace_back("armijo_passed",armijo?1.:0.);
                audit.measurements.emplace_back("root_tolerance_passed",root?1.:0.);
                status=armijo?"MODEL_AGREEMENT_UNRESOLVED":"ARMIJO_FAILED";
                if(prediction.size())
                {
                    audit.arrays.emplace_back("raw_remainder",evaluation.residual-prediction);
                    const double decrease=base-predicted,scaledDecrease=scaledBase-scaledPredicted;
                    if(decrease>0.)
                    {
                        const double ratio=(base-actual)/decrease;
                        if(std::isfinite(ratio))audit.measurements.emplace_back("physical_decrease_ratio",ratio);
                        useful=resolved&&armijo&&std::isfinite(ratio)&&ratio>=.1&&
                            decrease>4.*fraction*physicalAllowance+64.*std::numeric_limits<double>::epsilon()*base&&
                            scaledDecrease>4.*fraction*scaledAllowance+64.*std::numeric_limits<double>::epsilon()*scaledBase;
                    }
                    if(useful||root)status="DECREASE_OBSERVED";
                }
                useful=useful||root;
            }
            audit.notes.emplace_back("decision",status);
            cb.ObserveActionAudit(audit);
            if(useful)break;
        }
    }

    SteadyKrylovActionAudit CorrectionBasis(const Vec &x,const Vec &f,const Mat &z,bool bounded)
    {
        Eigen::ColPivHouseholderQR<Mat> qr(z);
        qr.setThreshold(1e-11);
        const int n=static_cast<int>(z.rows()),m=static_cast<int>(qr.rank());
        Require(m>0,"frozen correction subspace is empty");
        const Mat q=qr.householderQ()*Mat::Identity(n,m);
        Mat response(n,m+1);
        response.col(0)=f;
        for(int j=0;j<m;++j)
            response.col(j+1)=Jacobian(x,f,q.col(j),1.,"orthonormal-basis",-1,j);
        auto data=AuditData(bounded?"bounded-current-basis":"frozen-size-basis",x,f,Vec::Zero(n),0,m);
        const double target=o.forcingInitial*Norm(f);
        data.measurements={{"physical_target",target},{"scaled_target",o.forcingInitial*sr.cwiseProduct(f).stableNorm()},
            {"forcing",o.forcingInitial},{"basis_rank",double(m)},{"krylov_columns",double(z.cols())},
            {"scaled_basis_orthogonality",(q.transpose()*q-Mat::Identity(m,m)).norm()}};
        data.arrays.emplace_back("scaled_correction_basis",q);
        data.arrays.emplace_back("raw_residual_and_actions",response);
        // Normalize each column before polarization. The multiplier is
        // recorded separately so it cannot change either forcing target.
        Vec normalization(m+1);
        Mat normalized(n,m+1);
        std::vector<std::array<double,6>> diagonals;
        for(int j=0;j<=m;++j)
        {
            normalization[j]=sr.cwiseProduct(response.col(j)).stableNorm();
            if(normalization[j]==0.)normalization[j]=1.;
            normalized.col(j)=response.col(j)/normalization[j];
            diagonals.push_back(cb.ResidualBlocks(normalized.col(j)));
        }
        std::array<Mat,6> grams;
        for(auto &gram:grams)gram=Mat::Zero(m+1,m+1);
        for(int i=0;i<=m;++i)
        {
            for(std::size_t k=0;k<grams.size();++k)
                grams[k](i,i)=diagonals[i][k]*diagonals[i][k];
            for(int j=0;j<i;++j)
            {
                Deadline();
                const auto plus=cb.ResidualBlocks(normalized.col(i)+normalized.col(j));
                const auto minus=cb.ResidualBlocks(normalized.col(i)-normalized.col(j));
                for(std::size_t k=0;k<grams.size();++k)
                    grams[k](i,j)=grams[k](j,i)=(plus[k]-minus[k])*(plus[k]+minus[k])/4.;
            }
        }
        data.arrays.emplace_back("response_column_normalization",normalization);
        for(std::size_t k=0;k<grams.size();++k)
            data.arrays.emplace_back("normalized_physical_gram_"+std::to_string(k),std::move(grams[k]));
        for(const auto &stencil:out.finiteDifferences)
            if(stencil.newton==newton)data.stencils.push_back(stencil);
        return data;
    }

    void AuditSpacings(const char *stage,const Vec &x,const Vec &f,const Vec &direction,
                       int cycle,int iteration,double forcing)
    {
        auto data=AuditData(stage,x,f,direction,cycle,iteration);
        data.measurements={{"forcing",forcing},{"current_fd",o.fdRelativeStep}};
        const auto first=out.finiteDifferences.size();
        const auto began=Clock::now();
        std::vector<std::pair<double,Vec>> samples;
        Vec reference, middle, previous;
        for(double multiplier:{1.,.5,.25,.125,.5})
        {
            const std::string label=std::to_string(samples.size());
            std::pair<double,Vec> sample;
            Vec action=Jacobian(x,f,direction,multiplier,"spacing_"+label,cycle,iteration,&sample);
            if(samples.size()==3)reference=action;
            if(samples.size()==1)middle=action;
            if(samples.size()==4)
                data.arrays.emplace_back("raw_repeated_action_difference",action-middle);
            else if(!samples.empty())
                data.arrays.emplace_back("raw_spacing_change_"+label,action-previous);
            previous=action;
            data.arrays.emplace_back("raw_model_residual_"+label,f+action);
            data.arrays.emplace_back("raw_action_"+label,std::move(action));
            data.arrays.emplace_back("raw_difference_"+label,sample.second);
            data.measurements.emplace_back("signed_parameter_"+label,sample.first);
            samples.push_back(std::move(sample));
        }
        for(std::size_t i=0;i<samples.size();++i)
            data.arrays.emplace_back("raw_taylor_"+std::to_string(i),
                                     samples[i].second-samples[i].first*reference);
        if(cb.ColoredAction&&factorSerial)
        {
            const Vec colored=cb.ColoredAction(sx.cwiseProduct(direction));
            data.arrays.emplace_back("raw_colored_action",colored);
            data.arrays.emplace_back("raw_colored_minus_middle",colored-middle);
        }
        if(cb.TraceDirection)
        {
            const Vec trace=cb.TraceDirection(sx.cwiseProduct(direction)).cwiseQuotient(sx);
            const Vec complement=direction-trace;
            if(trace.cwiseAbs().maxCoeff()>0.&&complement.cwiseAbs().maxCoeff()>0.)
            {
                const Vec a=Jacobian(x,f,trace,.5,"matching-k",cycle,iteration);
                const Vec b=Jacobian(x,f,complement,.5,"remaining-free",cycle,iteration);
                data.arrays.emplace_back("matching_k_direction",sx.cwiseProduct(trace));
                data.arrays.emplace_back("raw_matching_k_action",a);
                data.arrays.emplace_back("raw_complement_action",b);
                data.arrays.emplace_back("raw_additivity_gap",middle-a-b);
                if(cb.ColoredAction&&factorSerial)
                    data.arrays.emplace_back("raw_colored_matching_k",cb.ColoredAction(sx.cwiseProduct(trace)));
            }
        }
        data.stencils.assign(out.finiteDifferences.begin()+first,out.finiteDifferences.end());
        data.measurements.emplace_back("audit_seconds",Seconds(began));
        cb.ObserveActionAudit(data);
    }

    void AuditGap(const Vec &x,const Vec &f,const Vec &candidate,const Vec &freshAction,
                  const Vec &d0,const Vec &action0,const Vec &r0,const Mat &v,const Mat &z,
                  const Mat &w,const Mat &h,const Vec &y,const SteadyLinearRecord &check,double forcing)
    {
        actionAuditCaptured=true;
        const int k=int(y.size()),basis=k+(check.breakdown?0:1);
        const auto started=Clock::now();
        auto data=AuditData("failed-krylov",x,f,candidate,check.cycle,check.iteration);
        const Vec stored=r0-w.leftCols(k)*y, arnoldi=r0-v.leftCols(k+1)*(h*y);
        const Vec fresh=-sr.cwiseProduct(f+freshAction);
        const Mat orthogonality=v.leftCols(basis).transpose()*v.leftCols(basis)-Mat::Identity(basis,basis);
        double directionSum=d0.stableNorm(),actionSum=sr.cwiseProduct(action0).stableNorm();
        Vec replay=sx.cwiseProduct(d0);
        for(int j=0;j<k;++j)
        {
            directionSum+=std::abs(y[j])*z.col(j).stableNorm();
            actionSum+=std::abs(y[j])*w.col(j).stableNorm();
            replay+=y[j]*sx.cwiseProduct(z.col(j));
        }
        data.measurements={{"forcing",forcing},{"projected_scaled",check.projectedNorm},
            {"rank",double(check.leastSquaresRank)},{"basis_columns",double(basis)},
            {"physical_target",forcing*Norm(f)},{"scaled_target",forcing*sr.cwiseProduct(f).stableNorm()},
            {"arnoldi_relation_frobenius",(w.leftCols(k)-v.leftCols(k+1)*h).norm()},
            {"orthogonality_frobenius",orthogonality.norm()},
            {"least_squares_reconstruction",std::abs(arnoldi.stableNorm()-check.projectedNorm)},
            {"physical_direction_reconstruction",(replay-sx.cwiseProduct(candidate)).stableNorm()},
            {"scaled_direction_sum",directionSum},{"scaled_direction_norm",candidate.stableNorm()},
            {"scaled_action_sum",actionSum},{"scaled_summed_action_norm",(r0-stored+sr.cwiseProduct(action0)).stableNorm()}};
        data.arrays.emplace_back("scaled_seed_direction",d0);
        data.arrays.emplace_back("raw_seed_action",action0);
        data.arrays.emplace_back("scaled_arnoldi_basis",v.leftCols(k+1));
        data.arrays.emplace_back("scaled_preconditioned_directions",z.leftCols(k));
        data.arrays.emplace_back("scaled_stored_actions",w.leftCols(k));
        data.arrays.emplace_back("hessenberg",h);
        data.arrays.emplace_back("least_squares_coefficients",y);
        data.arrays.emplace_back("scaled_restart_rhs",r0);
        data.arrays.emplace_back("raw_fresh_action",freshAction);
        data.arrays.emplace_back("raw_stored_residual",stored.cwiseQuotient(sr));
        data.arrays.emplace_back("raw_arnoldi_residual",arnoldi.cwiseQuotient(sr));
        data.arrays.emplace_back("raw_fresh_residual",fresh.cwiseQuotient(sr));
        data.arrays.emplace_back("raw_additivity_gap",(stored-fresh).cwiseQuotient(sr));
        if(cb.FactoredAction)
        {
            const Vec z0=sx.cwiseProduct(z.col(0)),rhs=v.col(0).cwiseQuotient(sr);
            data.arrays.emplace_back("raw_first_backsolve_residual",cb.FactoredAction(z0)-rhs);
            data.arrays.emplace_back("raw_first_unshifted_factor_residual",cb.ColoredAction(z0)-rhs);
        }
        data.measurements.emplace_back("audit_seconds",Seconds(started));
        cb.ObserveActionAudit(data);
        AuditSpacings("failed-direction",x,f,candidate,check.cycle,check.iteration,forcing);
        if(cb.ColoredAction)
        {
            // Preserve the owner's current Krylov work and forcing ledger.
            // The control uses the same assembled factor and original matrix.
            Runner control(cb,sx,sr,o);
            control.start=start;control.newton=newton;
            const auto fixedStart=Clock::now();
            const auto fixed=control.Linear(x,f,Norm(f),forcing,30,nullptr,true);
            auto fixedData=AuditData("fixed-colored-control",x,f,fixed.correction,check.cycle,check.iteration);
            fixedData.measurements={{"forcing",forcing},{"matrix_passed",fixed.converged?1.:0.},
                {"products",double(control.out.statistics.linearIterations)},{"seconds",Seconds(fixedStart)}};
            const auto last=std::find_if(control.out.linear.rbegin(),control.out.linear.rend(),
                                         [](const auto &line){return line.trueChecked;});
            if(last!=control.out.linear.rend())
                fixedData.measurements.emplace_back("projected_scaled_at_last_check",last->projectedNorm);
            fixedData.arrays.emplace_back("raw_colored_residual",f+fixed.action);
            cb.ObserveActionAudit(fixedData);
            if(fixed.correction.cwiseAbs().maxCoeff()>0.)
                AuditSpacings("fixed-colored-direction",x,f,fixed.correction,check.cycle,check.iteration,forcing);
        }
    }

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
            if(o.merit->bounded)
            {
                state.radius.emplace();
                state.radius->next=o.merit->bounded->initialRadius;
            }
        }
        out.meritStart=out.meritEnd=state;
        if(cb.ObserveMeritState)cb.ObserveMeritState(state);
    }

    void BoundedDirections(Vec &x,Vec &f,double &norm,int &factorBase)
    {
        const auto &policy=*o.merit->bounded;
        const auto starts=out.statistics.jacobianActions;
        if(out.meritEnd->radius->stagnated)
            throw Stop("BOUNDED_STAGNATION: reanchor on new numerical evidence");
        if(factorBase>=0&&newton-factorBase>o.maximumInverseAge)
            throw Stop("BOUNDED_INVERSE_AGE: measured run ended at its declared inverse age");
        if(cb.SetCorrectionMass)cb.SetCorrectionMass(0.);
        if(policy.recycledBasis.size())
        {
            boundedBasis=CorrectionBasis(x,f,policy.recycledBasis,true);
            cb.ObserveActionAudit(*boundedBasis);
        }
        else
            Correct(x,f,norm,o.forcingInitial,factorBase,factorBase<0?"bounded cold basis":"bounded current basis",
                    factorBase>=0,policy.columns);
        if(!boundedBasis)throw Stop("BOUNDED_BASIS: no current correction subspace");
        const auto q=std::find_if(boundedBasis->arrays.begin(),boundedBasis->arrays.end(),
                                  [](const auto &a){return a.first=="scaled_correction_basis";});
        const double scaledBase=sr.cwiseProduct(f).stableNorm();
        const double physicalTarget=o.forcingInitial*norm,scaledTarget=o.forcingInitial*scaledBase;
        const double roundoff=64.*std::numeric_limits<double>::epsilon();
        double radius=out.meritEnd->radius->next;
        if(q==boundedBasis->arrays.end())throw std::logic_error("bounded basis has no scaled coordinates");

        for(int proposal=0;proposal<policy.maximumProposals;++proposal)
        {
            Deadline();
            if(radius<policy.minimumRadius)break;
            const auto began=Clock::now();
            const auto proposalData=cb.OptimizeBounded(radius,proposal,o.maximumSeconds-Seconds(start));
            const double optimizerSeconds=Seconds(began);
            CheckVector(proposalData.coordinates,q->second.cols(),"bounded optimizer coordinates");
            const Vec direction=q->second*proposalData.coordinates;
            const double size=direction.cwiseAbs().maxCoeff();
            if(!(size<=radius*(1.+1e-10)))
                throw Stop("BOUNDED_RADIUS: optimizer result exceeds its declared radius");
            auto audit=AuditData("bounded-direction",x,f,direction,-1,0);
            audit.proposal=proposal;
            audit.notes.emplace_back("optimizer_status",proposalData.status);
            audit.arrays.emplace_back("optimizer_coordinates",proposalData.coordinates);
            audit.measurements={{"radius",radius},{"scaled_size",size},{"physical_base",norm},
                {"scaled_base",scaledBase},{"physical_target",physicalTarget},{"scaled_target",scaledTarget},
                {"optimizer_score",proposalData.score},{"optimizer_iterations",double(proposalData.iterations)},
                {"optimizer_seconds",optimizerSeconds}};
            if(proposalData.lower)
            {
                audit.measurements.emplace_back("optimizer_lower_bound",*proposalData.lower);
                audit.measurements.emplace_back("optimizer_gap",proposalData.score-*proposalData.lower);
            }
            const auto first=out.finiteDifferences.size();
            Vec primary,half,wide;
            double parameter=0.,halfParameter=0.;
            auto action=[&](double multiplier,const char *label,Vec &result,double *step=nullptr) {
                try
                {
                    result=Jacobian(x,f,direction,multiplier,label,-1,proposal);
                    if(step)*step=out.finiteDifferences.back().parameterStep;
                }
                catch(const Stop &error)
                {
                    if(!FDStop(error))throw;
                    audit.notes.emplace_back(label,error.what());
                }
            };
            std::string decision="NO_NONZERO_BOUNDED_PROPOSAL";
            if(size>0.)
            {
                action(1.,"bounded-primary",primary,&parameter);
                action(.5,"bounded-half",half,&halfParameter);
                action(2.,"bounded-wide",wide);
                decision="UNRESOLVED_CURRENT_ACTION";
            }
            Evaluation evaluation;
            Vec candidate;
            double trialNorm=0.,displacement=0.,ratio=0.;
            bool accepted=false;
            if(primary.size())
            {
                const Vec predicted=f+primary;
                const double model=Norm(predicted),scaledModel=sr.cwiseProduct(predicted).stableNorm();
                const double decrease=norm-model,scaledDecrease=scaledBase-scaledModel;
                audit.arrays.emplace_back("raw_action_primary",primary);
                audit.arrays.emplace_back("raw_model_residual",predicted);
                audit.measurements.emplace_back("physical_model",model);
                audit.measurements.emplace_back("scaled_model",scaledModel);
                audit.measurements.emplace_back("effective_physical_error_ratio",model/norm);
                audit.measurements.emplace_back("effective_scaled_error_ratio",scaledModel/scaledBase);
                audit.measurements.emplace_back("old_newton_model_pass",model<=physicalTarget&&scaledModel<=scaledTarget?1.:0.);
                if(wide.size())audit.arrays.emplace_back("raw_action_wide",wide);
                if(half.size()&&halfParameter<.8*parameter)
                {
                    const Vec allowance=2.*(half-primary);
                    const double physicalAllowance=Norm(allowance),scaledAllowance=sr.cwiseProduct(allowance).stableNorm();
                    audit.arrays.emplace_back("raw_action_half",half);
                    audit.arrays.emplace_back("raw_action_allowance",allowance);
                    audit.measurements.emplace_back("physical_action_allowance",physicalAllowance);
                    audit.measurements.emplace_back("scaled_action_allowance",scaledAllowance);
                    const bool resolved=physicalAllowance<.5*physicalTarget&&scaledAllowance<.5*scaledTarget;
                    const bool eligible=resolved&&decrease>4.*physicalAllowance+roundoff*norm&&
                        scaledDecrease>4.*scaledAllowance+roundoff*scaledBase;
                    decision=resolved?"INSUFFICIENT_RESOLVED_MODEL_DESCENT":decision;
                    if(eligible)
                    {
                        audit.notes.emplace_back("model_status",decrease>=.01*norm?"USEFUL_MODEL_DESCENT":"WEAK_MODEL_DESCENT");
                        candidate=x+sx.cwiseProduct(direction);
                        const Vec delta=(candidate-x).cwiseQuotient(sx);
                        displacement=delta.cwiseAbs().maxCoeff();
                        const double coordinateNorm=std::max(1.,x.cwiseQuotient(sx).cwiseAbs().maxCoeff());
                        if(displacement>radius*(1.+1e-10)+roundoff*coordinateNorm)
                            throw Stop("BOUNDED_RADIUS: represented trial exceeds its declared radius");
                        audit.arrays.emplace_back("trial_free",candidate);
                        audit.measurements.emplace_back("fraction",1.);
                        audit.measurements.emplace_back("scaled_displacement",displacement);
                        if(cb.ScaledUnknownBlocks)
                        {
                            const auto blocks=cb.ScaledUnknownBlocks(delta);
                            audit.arrays.emplace_back("scaled_displacement_blocks",Eigen::Map<const Vec>(blocks.data(),blocks.size()));
                        }
                        evaluation=Trial(candidate,false);
                        audit.measurements.emplace_back("nonlinear_validated",evaluation.nonlinearValidated?1.:0.);
                        decision="FINITE_STATE_REJECTED";
                        if(!evaluation.admissible)audit.notes.emplace_back("failure",evaluation.reason);
                        if(evaluation.admissible)
                        {
                            trialNorm=Norm(evaluation.residual);
                            const bool armijo=trialNorm<=(1.-o.armijo)*norm;
                            const bool root=trialNorm<=o.acceptanceNorm;
                            ratio=(norm-trialNorm)/decrease;
                            const bool agreement=std::isfinite(ratio)&&ratio>=.1;
                            accepted=root||(armijo&&agreement);
                            decision=accepted?"BOUNDED_ACCEPTED":(armijo?"INSUFFICIENT_ACTUAL_AGREEMENT":"ARMIJO_FAILED");
                            audit.arrays.emplace_back("raw_actual_residual",evaluation.residual);
                            audit.arrays.emplace_back("raw_remainder",evaluation.residual-predicted);
                            audit.measurements.emplace_back("physical_actual",trialNorm);
                            audit.measurements.emplace_back("scaled_actual",sr.cwiseProduct(evaluation.residual).stableNorm());
                            audit.measurements.emplace_back("armijo_limit",(1.-o.armijo)*norm);
                            audit.measurements.emplace_back("armijo_passed",armijo?1.:0.);
                            audit.measurements.emplace_back("ratio_passed",agreement?1.:0.);
                            audit.measurements.emplace_back("root_tolerance_passed",root?1.:0.);
                            if(std::isfinite(ratio))audit.measurements.emplace_back("physical_decrease_ratio",ratio);
                            Observe(candidate,evaluation.residual,trialNorm,false);
                        }
                    }
                }
            }
            audit.stencils.assign(out.finiteDifferences.begin()+first,out.finiteDifferences.end());
            const auto work=out.statistics.jacobianActions-starts;
            audit.measurements.emplace_back("current_products_at_base",double(work));
            audit.notes.emplace_back("decision",decision);
            cb.ObserveActionAudit(audit);
            Deadline();
            if(accepted)
            {
                AdvanceRadius(norm,trialNorm,radius,displacement,ratio,work,candidate-x);
                x=std::move(candidate);f=std::move(evaluation.residual);norm=trialNorm;
                Observe(x,f,norm,true);
                ++out.statistics.acceptedNewtonSteps;
                auto &record=out.newton.back();
                record.accepted=true;record.acceptedNorm=norm;record.acceptedStep=1.;
                boundedBasis.reset();
                if(out.meritEnd->radius->stagnated&&norm>o.acceptanceNorm)
                    throw Stop("BOUNDED_STAGNATION: six updates below one percent with no improving action rate");
                return;
            }
            radius*=.5;
        }
        throw Stop("BOUNDED_SEARCH_EXHAUSTED: no eligible actual original decrease in the radius budget");
    }

    void AdvanceRadius(double base,double actual,double used,double displacement,double ratio,
                       long long work,const Vec &delta)
    {
        auto state=*out.meritEnd;
        auto &radius=*state.radius;
        state.previousNorm=base;state.acceptedNorm=actual;state.step=1.;
        state.ratio=std::isfinite(ratio)?std::optional<double>(ratio):std::nullopt;
        ++state.accepted;
        radius.used=used;radius.displacement=displacement;radius.acceptedDisplacement=delta;
        radius.next=used;
        state.decision="hold radius";
        if(ratio>.75&&displacement>=.8*used)
        {
            radius.next=std::min(2.*used,o.merit->bounded->maximumRadius);
            if(radius.next>used)state.decision="grow radius";
        }
        const double rate=(base-actual)/base/std::max(1LL,work);
        if(!radius.windowUpdates){radius.windowBase=base;radius.firstRate=rate;}
        radius.lastRate=rate;
        if(++radius.windowUpdates==6)
        {
            radius.stagnated=(radius.windowBase-actual)/radius.windowBase<.01&&radius.lastRate<=radius.firstRate;
            if(!radius.stagnated)radius.windowUpdates=0;
        }
        out.meritEnd=state;
        out.meritAcceptances.push_back(state);
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

    void CaptureMerit(const Vec &x,const Vec &f,const Vec &direction,
                      const Vec &primary,const Vec &finer,const Vec &original,
                      const FDSample &primarySample,const FDSample &halfSample,
                      std::size_t firstStencil,const SteadyMeritRecord &rec)
    {
        if(!cb.ObserveActionAudit||!rec.linearPassed)return;
        const bool unresolved=rec.decision=="ACTION_UNRESOLVED";
        // A prospective first usable direction and up to two failures per solve.
        if(unresolved?meritFailureCaptures>=2:meritPanelCaptured)return;
        if(unresolved)++meritFailureCaptures;
        else meritPanelCaptured=true;
        auto data=AuditData(unresolved?"merit-unresolved":"merit-current-panel",x,f,direction,-1,-1);
        data.proposal=rec.proposal;
        data.notes={{"decision",rec.decision}};
        data.measurements={{"forcing",rec.forcing},{"physical_target",rec.forcing*rec.baseNorm},
            {"scaled_target",rec.forcing*rec.scaledBaseNorm},{"current_fd",o.fdRelativeStep},
            {"primary_parameter",rec.primaryParameter},{"half_parameter",rec.halfParameter}};
        if(rec.modelEvaluated)
        {
            data.measurements.emplace_back("physical_spacing_allowance",rec.actionUncertainty);
            data.measurements.emplace_back("scaled_spacing_allowance",rec.scaledUncertainty);
            data.measurements.emplace_back("active_slope",rec.activeSlope);
        }
        data.arrays.emplace_back("scaled_direction",direction);
        data.arrays.emplace_back("raw_primary_action",primary);
        data.arrays.emplace_back("raw_original_prediction",f+original);
        if(finer.size())data.arrays.emplace_back("raw_half_action",finer);

        // Same direction, factor and scales. Retain the actual production and
        // half-spacing samples, and repeat the production action independently.
        struct PanelSample { Vec action; FDSample sample; };
        std::array<PanelSample,5> panel;
        panel[1]={primary,primarySample};
        panel[2]={finer,halfSample};
        const std::array<double,5> multipliers{2.,1.,.5,.25,1.};
        for(int i:{0,3,4})
        {
            try
            {
                panel[i].action=Jacobian(x,f,direction,multipliers[i],
                    "merit-panel-"+std::to_string(i),-1,rec.proposal,&panel[i].sample);
            }
            catch(const Stop &error)
            {
                if(!FDStop(error))throw;
                data.notes.emplace_back("unavailable_sample_"+std::to_string(i),error.what());
            }
        }
        const Vec &reference=panel[3].action.size()?panel[3].action:primary;
        for(std::size_t i=0;i<panel.size();++i)
        {
            if(!panel[i].action.size())continue;
            const auto label=std::to_string(i);
            const auto &[h,difference]=panel[i].sample;
            data.measurements.emplace_back("signed_parameter_"+label,h);
            data.arrays.emplace_back("raw_action_"+label,panel[i].action);
            data.arrays.emplace_back("raw_difference_"+label,difference);
            data.arrays.emplace_back("raw_taylor_"+label,difference-h*reference);
            if(i>0&&i<4&&panel[i-1].action.size())
                data.arrays.emplace_back("raw_spacing_change_"+label,panel[i].action-panel[i-1].action);
        }
        if(panel[4].action.size())
            data.arrays.emplace_back("raw_repeated_action_difference",panel[4].action-primary);
        data.stencils.assign(out.finiteDifferences.begin()+firstStencil,out.finiteDifferences.end());
        cb.ObserveActionAudit(data);
    }

    void CaptureRejectedTrial(const Vec &x,const Vec &f,const Vec &d,const Vec &action,
                              const Vec &candidate,const Evaluation &evaluation,
                              const SteadyLineSearchRecord &trial)
    {
        if(!cb.ObserveActionAudit||rejectedTrialCaptured)return;
        rejectedTrialCaptured=true;
        const Vec prediction=f+trial.step*action;
        auto data=AuditData("merit-armijo-rejected",x,f,d.cwiseQuotient(sx),-1,trial.trial);
        data.proposal=trial.proposal;
        data.notes={{"decision","FEASIBLE_ARMIJO_REJECTED"}};
        data.measurements={{"fraction",trial.step},{"actual_norm",trial.trueNorm},
            {"predicted_norm",*trial.predictedNorm},{"scaled_displacement",trial.scaledDisplacement}};
        data.arrays.emplace_back("trial_free",candidate);
        data.arrays.emplace_back("raw_original_action",action);
        data.arrays.emplace_back("raw_plain_residual",evaluation.residual);
        data.arrays.emplace_back("raw_original_prediction",prediction);
        data.arrays.emplace_back("raw_coupled_remainder",evaluation.residual-prediction);
        cb.ObserveActionAudit(data);
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
        rec.scaledDirectionNorm=linear.correction.cwiseAbs().maxCoeff();
        if(cb.ScaledUnknownBlocks)rec.directionBlocks=cb.ScaledUnknownBlocks(linear.correction);
        const bool nonzero=rec.scaledDirectionNorm!=0.;
        const auto firstStencil=out.finiteDifferences.size();
        FDSample primarySample,halfSample;
        const Vec primary=nonzero?Jacobian(x,f,linear.correction,1.,"original-merit-primary",-1,proposal,
                                           cb.ObserveActionAudit?&primarySample:nullptr):
                                  linear.action;
        rec.primaryParameter=nonzero?out.finiteDifferences.back().parameterStep:0.;
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
        Vec finer;
        if(multiplier<1.)
        {
            try
            {
                finer=Jacobian(x,f,linear.correction,multiplier,"original-merit-spacing",-1,proposal,
                               cb.ObserveActionAudit?&halfSample:nullptr);
                rec.halfParameter=out.finiteDifferences.back().parameterStep;
            }
            catch(const Stop &error)
            {
                if(!FDStop(error)||!cb.ObserveActionAudit)throw;
            }
        }
        auto publish=[&] {
            CaptureMerit(x,f,linear.correction,primary,finer,original,
                         primarySample,halfSample,firstStencil,rec);
        };
        if(!finer.size()||rec.halfParameter>=.8*rec.primaryParameter)
        {
            publish();return rec;
        }
        const Vec change=finer-primary;
        rec.actionUncertainty=2.*Norm(change);
        rec.scaledUncertainty=2.*sr.cwiseProduct(change).stableNorm();
        rec.uncertaintyBlocks=cb.ResidualBlocks(2.*change);
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
        publish();
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
                double step=assessment.step;
                if(o.merit->initialTrialFraction>0.)
                {
                    const auto &previous=*out.meritEnd;
                    double start=o.merit->initialTrialFraction;
                    if(previous.step>0.)
                    {
                        const bool reliable=previous.ratio&&*previous.ratio>=.5&&*previous.ratio<=1.5;
                        start=(reliable?2.:1.)*previous.step;
                    }
                    step=std::min(step,start);
                }
                if(Backtrack(x,f,norm,PhysicalCorrection(result),original,step,
                             assessment.activeSlope,std::min(8,o.maximumLineSearchTrials),proposal,
                             &assessment,remaining()))
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

    std::optional<CurvedCandidate> CorrectTrial(
        const Vec &x,const Vec &f,double norm,const Vec &d,const Vec &action,
        const Vec &plainX,const Evaluation &plain,const SteadyLineSearchRecord &trial,
        const SteadyMeritRecord &primary,int &availableProducts)
    {
        const auto &policy=*o.merit->curvature;
        auto &rec=out.curvatureTrials.emplace_back();
        Timer timer(rec.seconds);
        rec.base=newton;rec.proposal=trial.proposal;rec.trial=trial.trial;rec.step=trial.step;
        rec.factorSerial=factorSerial;
        rec.plainNorm=trial.trueNorm;rec.primaryPrediction=*trial.predictedNorm;
        rec.armijoLimit=(1.-o.armijo*trial.step)*norm;
        const Vec predicted=f+trial.step*action;
        const Vec defect=plain.residual-predicted;
        rec.defectBlocks=*trial.remainderBlocks;
        const auto &primaryError=*primary.uncertaintyBlocks;
        const double roundoff=64.*std::numeric_limits<double>::epsilon()*norm;
        bool resolved=true;
        for(std::size_t i=0;i<rec.defectBlocks.size();++i)
        {
            rec.primaryUncertainty[i]=trial.step*primaryError[i];
            rec.blockBudgets[i]=policy.budgetFraction*
                (rec.armijoLimit-(*trial.predictedBlocks)[i]-rec.primaryUncertainty[i]-roundoff);
            if((*trial.actualBlocks)[i]>rec.armijoLimit&&
               rec.defectBlocks[i]<=2.*rec.primaryUncertainty[i]+roundoff)
                resolved=false;
        }
        const double scaledDefect=sr.cwiseProduct(defect).stableNorm();
        resolved=resolved&&scaledDefect>2.*trial.step*primary.scaledUncertainty;
        const double scaledSlack=trial.step*(primary.scaledBaseNorm-primary.scaledError-primary.scaledUncertainty);
        rec.scaledBudget=std::min(.25*scaledDefect,policy.budgetFraction*scaledSlack);
        rec.scaledSizeLimit=policy.maximumScaledSize*trial.scaledDisplacement;
        if(!resolved)
        {
            rec.decision="DEFECT_BELOW_PRIMARY_ACTION_RESOLUTION";
            return std::nullopt;
        }
        if(rec.scaledBudget<=0.||*std::min_element(rec.blockBudgets.begin(),rec.blockBudgets.end())<=0.)
        {
            rec.decision="NO_VERIFIED_MODEL_MARGIN";
            return std::nullopt;
        }
        const int cap=std::min(availableProducts,policy.maximumProducts);
        if(cap<4)
        {
            rec.decision="SHARED_WORK_CAP";
            return std::nullopt;
        }

        std::optional<SteadyKrylovActionAudit> receipt;
        const auto firstStencil=out.finiteDifferences.size();
        if(cb.ObserveActionAudit&&curvatureCaptures<2)
        {
            ++curvatureCaptures;
            receipt=AuditData("merit-curvature-trial",x,f,d.cwiseQuotient(sx),-1,trial.trial);
            receipt->proposal=trial.proposal;
            receipt->arrays.emplace_back("trial_free",plainX);
            receipt->arrays.emplace_back("raw_original_prediction",predicted);
            receipt->arrays.emplace_back("raw_plain_residual",plain.residual);
            receipt->arrays.emplace_back("raw_coupled_remainder",defect);
        }
        const auto productsBefore=out.statistics.jacobianActions;
        auto finish=[&](const char *decision,std::optional<CurvedCandidate> candidate=std::nullopt) {
            rec.decision=decision;
            rec.products=int(out.statistics.jacobianActions-productsBefore);
            availableProducts-=rec.products;
            if(receipt)
            {
                receipt->notes={{"decision",rec.decision},{"reason",rec.reason}};
                receipt->measurements={{"fraction",trial.step},{"plain_norm",rec.plainNorm},
                    {"armijo_limit",rec.armijoLimit},{"physical_target",*std::min_element(
                        rec.blockBudgets.begin(),rec.blockBudgets.end())},
                    {"scaled_target",rec.scaledBudget},{"products",double(rec.products)},
                    {"scaled_size_limit",rec.scaledSizeLimit}};
                receipt->stencils.assign(out.finiteDifferences.begin()+firstStencil,out.finiteDifferences.end());
                cb.ObserveActionAudit(*receipt);
            }
            return candidate;
        };
        // The RHS defect and F0(x) are distinct. In a one-sided stencil the
        // action still subtracts F0(x), never this auxiliary load.
        const AuxiliaryRhs rhs{defect,rec.blockBudgets,rec.scaledBudget,cap-2};
        LinearResult auxiliary;
        try
        {
            auxiliary=Linear(x,f,Norm(defect),.25,cap-2,nullptr,false,&rhs);
        }
        catch(const Stop &error)
        {
            if(!FDStop(error))throw;
            rec.reason=error.what();
            return finish("AUXILIARY_ACTION_UNAVAILABLE");
        }
        rec.scaledCorrection=auxiliary.correction.cwiseAbs().maxCoeff();
        if(*rec.scaledCorrection==0.)
            return finish("AUXILIARY_LINEAR_TARGET_NOT_MET");
        Vec auxiliaryAction,halfAction;
        const double multiplier=std::max(.5,o.fdMinimumRelativeStep/o.fdRelativeStep);
        if(multiplier>=1.)return finish("AUXILIARY_ACTION_UNRESOLVED");
        try
        {
            auxiliaryAction=Jacobian(x,f,auxiliary.correction,1.,"curvature-primary",-1,trial.trial);
            const double parameter=out.finiteDifferences.back().parameterStep;
            halfAction=Jacobian(x,f,auxiliary.correction,multiplier,"curvature-half",-1,trial.trial);
            if(out.finiteDifferences.back().parameterStep>=.8*parameter)
                return finish("AUXILIARY_ACTION_UNRESOLVED");
        }
        catch(const Stop &error)
        {
            if(!FDStop(error))throw;
            rec.reason=error.what();
            return finish("AUXILIARY_ACTION_UNAVAILABLE");
        }
        const Vec error=defect+auxiliaryAction, uncertainty=2.*(auxiliaryAction-halfAction);
        rec.errorBlocks=cb.ResidualBlocks(error);
        rec.uncertaintyBlocks=cb.ResidualBlocks(uncertainty);
        rec.scaledError=sr.cwiseProduct(error).stableNorm();
        rec.scaledUncertainty=sr.cwiseProduct(uncertainty).stableNorm();
        bool resolvedAuxiliary=*rec.scaledError+*rec.scaledUncertainty<=rec.scaledBudget;
        for(std::size_t i=0;i<rec.defectBlocks.size();++i)
            resolvedAuxiliary=resolvedAuxiliary&&
                (*rec.errorBlocks)[i]+(*rec.uncertaintyBlocks)[i]<=rec.blockBudgets[i];
        if(receipt)
        {
            receipt->arrays.emplace_back("auxiliary_scaled_correction",auxiliary.correction);
            receipt->arrays.emplace_back("raw_auxiliary_action",auxiliaryAction);
            receipt->arrays.emplace_back("raw_auxiliary_half_action",halfAction);
            receipt->arrays.emplace_back("raw_auxiliary_residual",error);
        }
        if(!auxiliary.converged)return finish("AUXILIARY_LINEAR_TARGET_NOT_MET");
        if(!resolvedAuxiliary)return finish("AUXILIARY_ACTION_AND_ERROR_OVER_BUDGET");
        if(*rec.scaledCorrection>rec.scaledSizeLimit)return finish("AUXILIARY_SIZE_LIMIT");
        Vec corrected=plainX+sx.cwiseProduct(auxiliary.correction);
        const Vec displacement=(corrected-x).cwiseQuotient(sx);
        rec.correctedDisplacement=displacement.cwiseAbs().maxCoeff();
        auto evaluation=Trial(corrected,false);
        rec.reason=evaluation.reason;
        if(!evaluation.admissible)return finish("CORRECTED_TRIAL_INADMISSIBLE");
        rec.correctedNorm=Norm(evaluation.residual);
        rec.correctedBlocks=cb.ResidualBlocks(evaluation.residual);
        rec.correctedPrediction=Norm(plain.residual+auxiliaryAction);
        const Vec remaining=evaluation.residual-predicted;
        rec.remainingDefectBlocks=cb.ResidualBlocks(remaining);
        if(receipt)
        {
            receipt->arrays.emplace_back("corrected_free",corrected);
            receipt->arrays.emplace_back("raw_corrected_residual",evaluation.residual);
            receipt->arrays.emplace_back("raw_remaining_defect",remaining);
        }
        if(*rec.correctedNorm>std::max(rec.armijoLimit,o.acceptanceNorm))
            return finish("CORRECTED_TRIAL_ARMIJO_REJECTED");
        return finish("CORRECTED_TRIAL_PASSED",CurvedCandidate{std::move(corrected),std::move(evaluation)});
    }

    bool Backtrack(Vec &x,Vec &f,double &norm,const Vec &correction,const Vec &originalAction,
                   double step,double modelSlope,int trialLimit,int proposal=-1,
                   const SteadyMeritRecord *assessment=nullptr,int availableProducts=0)
    {
        if((o.measureDirection||o.merit)&&cb.InspectDirection)cb.InspectDirection(newton,x,correction);
        bool accepted=false;
        int curvatureTrials=0;
        {
            Timer timer(out.statistics.lineSearchSeconds);
            for(int trial=0;trial<trialLimit;++trial)
            {
                Deadline();
                SteadyLineSearchRecord ls;
                ls.newton=newton;ls.trial=trial;ls.step=step;ls.proposal=proposal;
                Vec candidate=x+step*correction;
                if(!candidate.allFinite())
                    throw Stop("NONFINITE_NEWTON_TRIAL: solver arithmetic failed; no clipping or hidden retry");
                if((candidate.array()==x.array()).all())
                    throw Stop("NEWTON_ROUNDOFF: line-search correction no longer changes the state");
                const Vec displacement=(candidate-x).cwiseQuotient(sx);
                ls.scaledDisplacement=displacement.cwiseAbs().maxCoeff();
                if(cb.ScaledUnknownBlocks)ls.displacementBlocks=cb.ScaledUnknownBlocks(displacement);
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
                        if(o.merit)CaptureRejectedTrial(x,f,correction,originalAction,candidate,evaluation,ls);
                        if(assessment&&o.merit->curvature&&curvatureTrials<2&&availableProducts>=4)
                        {
                            ++curvatureTrials;
                            auto curved=CorrectTrial(x,f,norm,correction,originalAction,candidate,
                                evaluation,ls,*assessment,availableProducts);
                            if(curved)
                            {
                                candidate=std::move(curved->x);evaluation=std::move(curved->evaluation);
                                const auto &record=out.curvatureTrials.back();
                                ls.trueNorm=*record.correctedNorm;
                                ls.actualBlocks=record.correctedBlocks;
                                ls.remainderBlocks=record.remainingDefectBlocks;
                                ls.scaledDisplacement=*record.correctedDisplacement;
                                if(cb.ScaledUnknownBlocks)
                                    ls.displacementBlocks=cb.ScaledUnknownBlocks((candidate-x).cwiseQuotient(sx));
                                ls.sufficientDecrease=true;ls.curvatureCorrected=true;
                                ls.reason="coupled trial passed original norm-Armijo";
                                const double decrease=norm-*ls.predictedNorm;
                                ls.decreaseRatio.reset();
                                const double ratio=decrease>0.?(norm-ls.trueNorm)/decrease:0.;
                                if(decrease>0.&&std::isfinite(ratio))ls.decreaseRatio=ratio;
                            }
                        }
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
                const auto began=Clock::now();
                Timer timer(out.statistics.preparationSeconds);
                cb.PreparePreconditioner(x,f);
                record.preparationSeconds=Seconds(began);
                factorBase=newton;
                record.prepared=true;
                factorSerial=out.statistics.preconditionerPreparations;
                factorMass=currentMass;
            }
            out.newton.back().missingLinearReason="current linear solve";
            record.serial=factorSerial;
            if(o.merit)record.mass=factorMass;
            const auto began=Clock::now();
            auto result=Linear(x,f,norm,forcing,workLimit>=0?workLimit:(reuse?o.reuseLinearIterations:0),guess);
            record.linearSeconds=Seconds(began);
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
                 double multiplier,const std::string &purpose,int cycle,int iteration,
                 std::pair<double,Vec> *sample=nullptr)
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
        auto difference=[&](const Vec &high,const Vec &low,const Vec &fh,const Vec &fl)->Vec {
            Vec value=fh-fl;
            if(cb.ExactPartResidual)
            {
                const Vec splitHigh=cb.ExactPartResidual(high),splitLow=cb.ExactPartResidual(low);
                CheckVector(splitHigh,sx.size(),"exact partial residual");
                CheckVector(splitLow,sx.size(),"exact partial residual");
                value-=splitHigh-splitLow;
            }
            return value;
        };
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
            record.plusDisplacement=(plus-x).cwiseQuotient(sx).cwiseAbs().maxCoeff();
            record.minusDisplacement=(minus-x).cwiseQuotient(sx).cwiseAbs().maxCoeff();
            out.finiteDifferences.push_back(record);
            auto &audit=out.finiteDifferences.back();
            auto positive=Trial(plus,true);
            audit.plusAdmissible=positive.admissible;audit.plusReason=positive.reason;
            auto negative=Trial(minus,true);
            audit.minusAdmissible=negative.admissible;audit.minusReason=negative.reason;
            Vec action;
            if(positive.admissible&&negative.admissible)
            {
                action=difference(plus,minus,positive.residual,negative.residual)/(2*h);
                audit.method="central";
            }
            else if(positive.admissible||negative.admissible)
            {
                const bool forward=positive.admissible;
                const Vec half=x+(forward?.5:-.5)*h*direction;
                if((half.array()==x.array()).all())
                    throw Stop("FD_ROUNDOFF: feasible one-sided half step is not representable");
                audit.halfDisplacement=(half-x).cwiseQuotient(sx).cwiseAbs().maxCoeff();
                auto middle=Trial(half,true);
                audit.halfAdmissible=middle.admissible;audit.halfReason=middle.reason;
                if(middle.admissible)
                {
                    // Difference from base first to avoid a large -3*f term.
                    if(forward)action=(4*difference(half,x,middle.residual,f)-
                                        difference(plus,x,positive.residual,f))/h;
                    else action=(4*difference(x,half,f,middle.residual)-
                                 difference(x,minus,f,negative.residual))/h;
                    audit.method=forward?"forward-second-order":"backward-second-order";
                }
            }
            audit.seconds=Seconds(attemptStart);
            if(action.size()!=0)
            {
                if(sample)
                {
                    sample->first=positive.admissible?h:-h;
                    sample->second=(positive.admissible?positive.residual:negative.residual)-f;
                }
                if(cb.ExactPartAction)
                {
                    const Vec tangent=cb.ExactPartAction(x,direction);
                    CheckVector(tangent,sx.size(),"exact partial derivative");
                    action+=tangent;
                }
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

    LinearResult Linear(const Vec &x,const Vec &f,double baseNorm,double forcing,int workLimit,
                        const Vec *guess,bool fixed=false,const AuxiliaryRhs *auxiliary=nullptr)
    {
        Timer timer(out.statistics.linearSeconds);
        const Vec b=-sr.cwiseProduct(auxiliary?auxiliary->defect:f);
        CheckVector(b,sx.size(),"scaled Newton RHS");
        const double bNorm=b.stableNorm();
        if(!Positive(bNorm))throw Stop("LINEAR_SCALING: nonzero true residual has zero/nonfinite scaled norm");
        const int n=static_cast<int>(sx.size()),m=std::min(n,o.restartDimension);
        Vec best=guess?*guess:Vec::Zero(n).eval(),bestAction=Vec::Zero(n);
        double bestNorm=baseNorm,bestScaled=bNorm;
        int work=0;
        const bool bounded=o.merit&&o.merit->bounded;
        const bool audit=!fixed&&!auxiliary&&!o.frozenSizeColumns&&!bounded&&
            bool(cb.ObserveActionAudit)&&!actionAuditCaptured;
        const double scaledTarget=auxiliary?auxiliary->scaledBudget:forcing*bNorm;
        const auto initialProducts=out.statistics.jacobianActions;
        auto canSpend=[&](int needed) {
            return !auxiliary||out.statistics.jacobianActions-initialProducts+needed<=auxiliary->maximumProducts;
        };
        auto passes=[&](const SteadyLinearRecord &check) {
            if(!(check.trueScaledNorm<=scaledTarget))return false;
            if(!auxiliary)return check.trueNorm<=forcing*baseNorm;
            for(std::size_t i=0;i<auxiliary->blockBudgets.size();++i)
                if(!((*check.trueBlocks)[i]<=auxiliary->blockBudgets[i]))return false;
            return true;
        };
        auto product=[&](const Vec &direction,const char *purpose,int cycle,int iteration)->Vec {
            if(!fixed)return Jacobian(x,f,direction,1.,purpose,cycle,iteration);
            const Vec action=cb.ColoredAction(sx.cwiseProduct(direction));
            CheckVector(action,n,"retained unshifted colored action");
            return action;
        };
        for(int cycle=0;cycle<o.maximumLinearCycles;++cycle)
        {
            Deadline();++out.statistics.linearCycles;
            const Vec baseCorrection=best;
            Vec baseAction=Vec::Zero(n);
            Vec r=b;
            if(baseCorrection.cwiseAbs().maxCoeff()!=0)
            {
                if(!canSpend(1))break;
                baseAction=product(baseCorrection,"restart-residual",cycle,-1);
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
            restart.auxiliary=bool(auxiliary);
            restart.trueScaledNorm=beta;restart.trueNorm=initialTrue;
            if(cb.ResidualBlocks)restart.trueBlocks=cb.ResidualBlocks(raw);
            restart.converged=passes(restart);
            out.linear.push_back(restart);
            if(restart.converged)
                return {baseCorrection,true,initialTrue,std::move(baseAction)};
            if(!Positive(beta))throw Stop("LINEAR_SCALING: inconsistent zero/nonfinite Krylov residual");
            Mat v=Mat::Zero(n,m+1),z=Mat::Zero(n,m),h=Mat::Zero(m+1,m);
            Mat stored;
            if(audit&&!actionAuditCaptured)stored.resize(n,m);
            v.col(0)=r/beta;
            const int available=workLimit?std::min(m,workLimit-work):m;
            for(int j=0;j<available;++j)
            {
                if(!canSpend(2))return {best,false,bestNorm,std::move(bestAction)};
                Deadline();++out.statistics.linearIterations;++work;
                z.col(j)=Precondition(v.col(j));
                Vec w=sr.cwiseProduct(product(z.col(j),"arnoldi",cycle,j));
                CheckVector(w,n,"scaled Arnoldi action");
                if(stored.size())stored.col(j)=w;
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
                if(o.frozenSizeColumns||bounded)
                {
                    if(k==available||breakdown)
                    {
                        auto basis=CorrectionBasis(x,f,z.leftCols(k),bounded);
                        cb.ObserveActionAudit(basis);
                        if(bounded)boundedBasis=std::move(basis);
                        return {best,false,bestNorm,std::move(bestAction)};
                    }
                    continue;
                }
                const Mat small=h.topLeftCorner(k+1,k);
                Vec target=Vec::Zero(k+1);target[0]=beta;
                // Rank-revealing small dense QR also handles happy/singular
                // Arnoldi breakdown without dividing by a zero diagonal.
                Eigen::CompleteOrthogonalDecomposition<Mat> qr(small);
                const Vec coefficients=qr.solve(target);
                CheckVector(coefficients,k,"small Krylov least-squares solution");
                SteadyLinearRecord record;record.newton=newton;record.cycle=cycle;
                record.auxiliary=bool(auxiliary);
                record.iteration=j;record.leastSquaresRank=static_cast<int>(qr.rank());
                record.projectedNorm=(target-small*coefficients).stableNorm();record.breakdown=breakdown;
                if(breakdown||k==available||k%o.trueLinearCheckInterval==0||
                   record.projectedNorm<=scaledTarget||!canSpend(3))
                {
                    const Vec candidate=baseCorrection+z.leftCols(k)*coefficients;
                    // This is a new nonlinear-residual FD evaluation in the
                    // accumulated direction, not the Arnoldi recurrence norm.
                    Vec action=product(candidate,"true-linear-residual",cycle,j);
                    const Vec check=b-action.cwiseProduct(sr);
                    CheckVector(check,n,"recomputed linear residual");
                    record.trueChecked=true;record.trueScaledNorm=check.stableNorm();
                    const Vec rawCheck=(check.array()/sr.array()).matrix();
                    record.trueNorm=Norm(rawCheck);
                    if(cb.ResidualBlocks)record.trueBlocks=cb.ResidualBlocks(rawCheck);
                    record.converged=passes(record);
                    if(record.trueNorm<bestNorm||
                       (record.trueNorm==bestNorm&&record.trueScaledNorm<bestScaled))
                    {best=candidate;bestAction=action;bestNorm=record.trueNorm;bestScaled=record.trueScaledNorm;}
                    out.linear.push_back(record);Deadline();
                    if(audit&&!actionAuditCaptured&&!record.converged&&
                       candidate.cwiseAbs().maxCoeff()>0.&&
                       (record.projectedNorm<=forcing*bNorm||k==available||breakdown))
                        AuditGap(x,f,candidate,action,baseCorrection,baseAction,r,v,z,stored,small,
                                 coefficients,record,forcing);
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
    bool actionAuditCaptured=false;
    bool meritPanelCaptured=false,rejectedTrialCaptured=false;
    int meritFailureCaptures=0;
    int curvatureCaptures=0;
    long long factorSerial=0;
    double currentMass=0.,factorMass=0.;
    std::optional<SteadyKrylovActionAudit> boundedBasis;
    SteadyKrylovResult out;
};
}

std::string SteadyKrylovActionAudit::FileSuffix(std::size_t sequence) const
{
    if(newton<0)return stage; // Standalone probes retain their one-shot names.
    return stage+"-base-"+std::to_string(newton)+"-proposal-"+std::to_string(proposal)+
           "-sequence-"+std::to_string(sequence);
}

std::string SteadyMeritState::ToJson()const
{
    namespace J=RANSWallSolveJson;
    std::vector<std::pair<std::string,std::string>> fields{{"schema",J::Quote(Identity)},{"used_mass",J::Number(usedMass)},
        {"next_mass",J::Number(nextMass)},{"reference_norm",J::Number(referenceNorm)},
        {"previous_norm",J::Number(previousNorm)},{"accepted_norm",J::Number(acceptedNorm)},
        {"used_forcing",J::Number(usedForcing)},{"next_forcing",J::Number(nextForcing)},
        {"step",J::Number(step)},{"actual_to_predicted_decrease",ratio?J::Number(*ratio):"null"},
        {"accepted",J::Number(accepted)},{"good_steps",J::Number(goodSteps)},
        {"decision",J::Quote(decision)}};
    if(radius)
    {
        const auto &r=*radius;
        fields.emplace_back("radius",J::Object({
            {"schema",J::Quote(SteadyBoundedOptions::Identity)},
            {"used",J::Number(r.used)},{"next",J::Number(r.next)},
            {"scaled_displacement",J::Number(r.displacement)},{"window_base",J::Number(r.windowBase)},
            {"window_updates",J::Number(r.windowUpdates)},{"first_rate",J::Number(r.firstRate)},
            {"last_rate",J::Number(r.lastRate)},{"stagnated",J::Boolean(r.stagnated)}}));
    }
    return J::Object(fields);
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
    Require(bool(m_callbacks.ExactPartAction)==bool(m_callbacks.ExactPartResidual),
            "exact partial differentiation needs matched value and tangent callbacks");
    Require(o.frozenSizeColumns>=0&&o.frozenSizeColumns<=o.restartDimension&&
            (!o.frozenSizeColumns||(bool(m_callbacks.ObserveActionAudit)&&bool(m_callbacks.ResidualBlocks))),
            "frozen size bound requires a scoped audit observer and block norms");
    Require(o.frozenFiniteTrialLimit>=0&&o.frozenFiniteTrialLimit<=4&&
            (!o.frozenFiniteTrialLimit||(o.evaluateInitialOnly&&bool(m_callbacks.ValidateNonlinearState)&&
              bool(m_callbacks.ObserveActionAudit)&&!m_callbacks.ObserveAccepted)),
            "frozen finite trials require evaluate-only validation and nonaccepted audit custody");
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
                m.maximumWork>0&&std::isfinite(m.initialTrialFraction)&&
                m.initialTrialFraction>=0.&&m.initialTrialFraction<=1.,"invalid original-merit correction policy");
        if(m.bounded)
        {
            const auto &b=*m.bounded;
            Require(!m.curvature&&m.initialMass==0.&&m.reanchorMass==0.&&m.initialTrialFraction==0.&&
                    !o.evaluateInitialOnly&&!o.frozenSizeColumns&&bool(m_callbacks.OptimizeBounded)&&
                    bool(m_callbacks.ObserveActionAudit)&&bool(m_callbacks.ValidateNonlinearState)&&
                    Positive(b.minimumRadius)&&Positive(b.initialRadius)&&Positive(b.maximumRadius)&&
                    b.minimumRadius<=b.initialRadius&&b.initialRadius<=b.maximumRadius&&
                    b.columns>=1&&b.columns<=24&&b.columns<=o.restartDimension&&
                    (o.maximumInverseAge==0||b.columns<=o.reuseLinearIterations)&&
                    b.maximumProposals>0&&b.maximumProposals<=8&&
                    2*b.columns+3*b.maximumProposals<=m.maximumWork,
                    "invalid original radius-bounded direction policy");
            if(b.recycledBasis.size())
                Require(b.recycledBasis.rows()==m_unknownScale.size()&&b.recycledBasis.cols()==b.columns&&
                        b.recycledBasis.allFinite()&&
                        (b.recycledBasis.transpose()*b.recycledBasis-Mat::Identity(b.columns,b.columns)).norm()<1e-10,
                        "recycled directions must form the signed orthonormal scaled basis");
        }
        if(m.curvature)
            Require(m.initialMass==0.&&m.reanchorMass==0.&&m.curvature->maximumProducts>=4&&
                    m.curvature->maximumProducts<=m.maximumWork&&
                    Fraction(m.curvature->maximumScaledSize)&&Fraction(m.curvature->budgetFraction)&&
                    (!m.restored||(m.restored->usedMass==0.&&m.restored->nextMass==0.)),
                    "coupled trial curvature correction requires an unshifted policy and strict bounded work");
        if(m.restored)
        {
            const auto &r=*m.restored;
            Require(std::isfinite(r.usedMass)&&r.usedMass>=0.&&std::isfinite(r.nextMass)&&r.nextMass>=0.&&
                    Positive(r.referenceNorm)&&Positive(r.previousNorm)&&Positive(r.acceptedNorm)&&
                    r.usedForcing>=o.forcingMinimum&&r.usedForcing<=o.forcingMaximum&&
                    r.nextForcing>=o.forcingMinimum&&r.nextForcing<=o.forcingMaximum&&
                    r.accepted>=0&&r.goodSteps>=0&&r.goodSteps<2&&std::isfinite(r.step)&&r.step>=0.&&r.step<=1.&&
                    (!r.ratio||std::isfinite(*r.ratio))&&!r.decision.empty(),"invalid restored original-merit state");
            Require(bool(r.radius)==bool(m.bounded),"restored radius policy differs");
            if(r.radius)
            {
                const auto &b=*m.bounded;
                const auto &a=*r.radius;
                Require(r.usedMass==0.&&r.nextMass==0.&&a.next>=b.minimumRadius&&a.next<=b.maximumRadius&&
                        a.windowUpdates>=0&&a.windowUpdates<=6&&std::isfinite(a.windowBase)&&
                        std::isfinite(a.firstRate)&&std::isfinite(a.lastRate)&&
                        (r.accepted>0?(a.used>=b.minimumRadius&&a.used<=b.maximumRadius&&Positive(a.displacement)&&
                            a.acceptedDisplacement.size()==m_unknownScale.size()&&a.acceptedDisplacement.allFinite()&&
                            std::abs(a.acceptedDisplacement.cwiseQuotient(m_unknownScale).cwiseAbs().maxCoeff()-
                                     a.displacement)<=1e-12*a.displacement):
                            (a.used==0.&&a.displacement==0.&&!a.acceptedDisplacement.size())),
                        "restored radius or accepted displacement differs");
            }
        }
    }

}

SteadyKrylovResult SteadyKrylov::Solve(const SteadyVector &initial)const
{return Runner(m_callbacks,m_unknownScale,m_residualRowScale,m_options).Solve(initial);}

SteadyDerivativeProbe SteadyKrylov::ProbeDirection(const SteadyVector &base,
    const SteadyVector &direction,const std::vector<double> &multipliers)const
{return Runner(m_callbacks,m_unknownScale,m_residualRowScale,m_options).Probe(base,direction,multipliers);}

void SteadyKrylov::AuditDirection(const SteadyVector &base,const SteadyVector &physicalDirection)const
{
    Require(bool(m_callbacks.ObserveActionAudit),"fixed-direction audit requires a publisher");
    Runner(m_callbacks,m_unknownScale,m_residualRowScale,m_options).Audit(base,physicalDirection);
}

void SteadyKrylov::AuditDirections(const SteadyVector &base,const Eigen::MatrixXd &physicalDirections)const
{
    Require(bool(m_callbacks.ObserveActionAudit)&&physicalDirections.rows()==base.size()&&
            physicalDirections.cols()>1&&physicalDirections.allFinite(),"invalid frozen candidate batch");
    Runner(m_callbacks,m_unknownScale,m_residualRowScale,m_options).AuditDirections(base,physicalDirections);
}

void SteadyKrylov::AuditColumn(const SteadyVector &base,const SteadyVector &physicalDirection,
                              const SteadyColoredColumnProbe &column,const std::string &label)const
{
    Require(bool(m_callbacks.ObserveActionAudit),"column audit requires a publisher");
    Runner(m_callbacks,m_unknownScale,m_residualRowScale,m_options).AuditColumn(base,physicalDirection,column,label);
}
}
