// SPDX-License-Identifier: MIT
#ifndef NEKTAR_WJ_NATIVE_STEADY_KRYLOV_H
#define NEKTAR_WJ_NATIVE_STEADY_KRYLOV_H

#include <Eigen/Core>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

namespace Nektar::WJNative
{
using SteadyVector=Eigen::VectorXd;

/** The only residual/feasibility exception interpreted as a rejected trial.
 * At the initial/base state this is fatal. Other callback exceptions always
 * propagate: scientific/input/programming failures are not line-search events.
 */
class SteadyInadmissibleTrial : public std::runtime_error
{
public:
    explicit SteadyInadmissibleTrial(const std::string &reason)
        :std::runtime_error(reason){}
};

struct SteadyTrialValidity
{
    bool admissible=true;
    std::string reason;
};

struct SteadyKrylovCallbacks
{
    // All callbacks use PHYSICAL coordinates/residuals/corrections. Residual
    // must be the caller's current fixed equations with essential Gamma
    // eliminated. Steady preparation uses alpha=0; an implicit physical step
    // supplies its exact BDF mass/history. This engine adds no mass or time.
    std::function<SteadyVector(const SteadyVector &)> Residual;
    // A fixed norm of raw residuals, including all required blocks. The caller
    // owns its physical thresholds; return a norm ratio if acceptanceNorm=1.
    // This must be a norm (including sign symmetry and homogeneity), not a
    // residual-dependent denominator or a cached residual evaluation.
    std::function<double(const SteadyVector &)> TrueNorm;
    // Approximate J^{-1} applied to a raw residual RHS, returning physical dx.
    std::function<SteadyVector(const SteadyVector &)> Precondition;
    std::function<SteadyTrialValidity(const SteadyVector &)> ValidTrial;
    // Optional explicit rebuild/refresh at each Newton base. Residual trials
    // must not silently redefine the prepared base operator or norm weights.
    std::function<void(const SteadyVector &,const SteadyVector &)> PreparePreconditioner;
    // Optional complete scientific guard for actual nonlinear states. Called
    // after ValidTrial and before Residual/Observe at every initial/re-evaluated
    // base and nonlinear line-search candidate. NEVER called on FD samples:
    // ValidTrial and Residual still enforce their constitutive/domain guards.
    // Throw SteadyInadmissibleTrial to reject a candidate (fatal at a base);
    // all unrelated exceptions propagate. No generic validation cache exists.
    // A caller cache must match the complete x exactly and must not be replaced
    // by a merely domain-admissible FD evaluation. Absence preserves V1 behavior.
    std::function<void(const SteadyVector &)> ValidateNonlinearState;
};

struct SteadyKrylovOptions
{
    // Required explicit run policy; zero defaults deliberately reject use.
    double maximumSeconds=0, acceptanceNorm=0;
    int maximumNewtonIterations=50;
    int restartDimension=40, maximumLinearCycles=8;
    int trueLinearCheckInterval=5, maximumLineSearchTrials=24;
    int maximumFDReductions=8;
    double fdRelativeStep=1e-6, fdMinimumRelativeStep=1e-12;
    double fdReduction=.5, lineSearchReduction=.5, armijo=1e-4;
    double forcingInitial=.3, forcingMinimum=1e-3, forcingMaximum=.5;
    double forcingGamma=.9, forcingExponent=1.5;
    double arnoldiBreakdownRelative=1e-14;
    // If explicitly enabled, an incomplete linear solve may enter line search
    // only after its recomputed true linear norm improves on the base norm.
    bool allowIncompleteLinearStep=false;
};

struct SteadyKrylovStatistics
{
    long long residualEvaluations=0, validityEvaluations=0, normEvaluations=0;
    long long preconditionerApplications=0, preconditionerPreparations=0;
    long long jacobianActions=0, fdRejectedTrials=0, lineSearchRejectedTrials=0;
    long long linearIterations=0, linearCycles=0, acceptedNewtonSteps=0;
    long long nonlinearValidationEvaluations=0, nonlinearValidationRejectedTrials=0;
    double residualSeconds=0, validitySeconds=0, normSeconds=0;
    double preconditionerSeconds=0, preparationSeconds=0, arnoldiSeconds=0;
    double nonlinearValidationSeconds=0;
    // Inclusive timings: linear/line-search/total must not be summed together
    // with callback timings to infer elapsed time or cost.
    double linearSeconds=0, lineSearchSeconds=0, totalSeconds=0;
};

struct SteadyFDRecord
{
    int newton=-1, cycle=-1, iteration=-1, reduction=0;
    std::string purpose, method, plusReason, minusReason, halfReason;
    double relativeStep=0, parameterStep=0, scaledCoordinateStep=0;
    double scaledDirectionInfinityNorm=0, seconds=0;
    // Admissible means ValidTrial + Residual domain checks only. These flags
    // do not certify optional full nonlinear-state/principal validation.
    bool plusAdmissible=false, minusAdmissible=false, halfAdmissible=false;
};

struct SteadyLinearRecord
{
    int newton=0, cycle=0, iteration=0, leastSquaresRank=0;
    double projectedNorm=0, trueScaledNorm=0, trueNorm=0;
    bool trueChecked=false, breakdown=false, converged=false;
};

struct SteadyLineSearchRecord
{
    int newton=0, trial=0;
    double step=0, trueNorm=0;
    bool admissible=false, evaluated=false, sufficientDecrease=false, accepted=false;
    bool nonlinearValidationAttempted=false, nonlinearValidated=false;
    std::string reason;
};

struct SteadyNewtonRecord
{
    int iteration=0;
    double baseNorm=0, forcing=0, linearTrueNorm=0, acceptedNorm=0, acceptedStep=0;
    bool linearConverged=false, accepted=false;
};

struct SteadyKrylovResult
{
    bool converged=false, hasEvaluatedState=false, bestFromAcceptedIterate=true;
    bool nonlinearValidationEnabled=false;
    std::string status="NOT_CONVERGED", reason;
    // Best evaluated admissible state, with optional full nonlinear validation
    // passed, not necessarily an Armijo-accepted
    // iterate on failure. The caller must explicitly install any returned x.
    // If the deadline expires before an initial norm can be evaluated,
    // hasEvaluatedState=false and trueNorm/residual are not evidence.
    SteadyVector x, residual;
    double trueNorm=0, lastAcceptedNorm=0;
    SteadyKrylovStatistics statistics;
    std::vector<SteadyFDRecord> finiteDifferences;
    std::vector<SteadyLinearRecord> linear;
    std::vector<SteadyLineSearchRecord> lineSearch;
    std::vector<SteadyNewtonRecord> newton;
};

struct SteadyDerivativeProbe
{
    bool completed=false;
    bool nonlinearValidationEnabled=false;
    std::string reason;
    // Raw physical J(x)*direction at fixed, prospectively supplied FD scales.
    // Defaults 1,1/2,1/4 provide two successive halvings. Comparisons are
    // diagnostic; this method does not certify derivative accuracy.
    std::vector<double> stepMultipliers;
    std::vector<SteadyVector> actions;
    std::vector<double> successiveRelativeScaledChanges;
    SteadyKrylovStatistics statistics;
    std::vector<SteadyFDRecord> finiteDifferences;
};

/** Steady feasible Newton + restarted flexible right-preconditioned GMRES.
 * Unknown scaling: x=diag(unknownScale)*y. Residual scaling:
 * r_scaled=diag(residualRowScale)*Residual(x). Both vectors are strictly
 * positive and fixed throughout Solve. No physical tolerance is inferred.
 *
 * FD is central whenever both trials are admissible. A second-order one-sided
 * formula or step reduction is allowed only after explicit callback rejection.
 * In scaled direction v, h=fdRelativeStep*max(1,|x/sx|_inf)/|v|_inf,
 * and physical trials are x +/- h*diag(sx)*v. No FD step is selected from a
 * successful/failing residual magnitude. ProbeDirection's independent fixed
 * multipliers are for prospective Linux accuracy checks, not hidden adaptation.
 * Arnoldi uses two MGS passes. A newly evaluated FD action of the accumulated
 * correction checks the linear residual, independently of its Hessenberg
 * estimate. Newton accepts the current steady TrueNorm by feasible norm-Armijo
 * decrease, or direct satisfaction of acceptanceNorm. There is no clipping.
 * Optional full nonlinear validation precedes any candidate observation or
 * acceptance. FD points are uncommitted domain evaluations, never best states.
 *
 * The finite wall deadline is cooperative between operations/callbacks; it
 * cannot interrupt a blocking user callback. All algorithm cycle bounds are
 * declared options, not experimental-attempt limits. Base/input errors throw;
 * bounded algorithm failure/deadline returns NOT_CONVERGED with evidence.
 */
class SteadyKrylov
{
public:
    static constexpr const char *Identity="WJH27SteadyJFNKV1";
    SteadyKrylov(SteadyKrylovCallbacks callbacks,
                 SteadyVector unknownScale,SteadyVector residualRowScale,
                 SteadyKrylovOptions options);
    const SteadyKrylovOptions &Options()const{return m_options;}
    SteadyKrylovResult Solve(const SteadyVector &initial)const;
    SteadyDerivativeProbe ProbeDirection(const SteadyVector &base,
        const SteadyVector &physicalDirection,
        const std::vector<double> &stepMultipliers={1.,.5,.25})const;
private:
    SteadyKrylovCallbacks m_callbacks;
    SteadyVector m_unknownScale,m_residualRowScale;
    SteadyKrylovOptions m_options;
};
}
#endif
