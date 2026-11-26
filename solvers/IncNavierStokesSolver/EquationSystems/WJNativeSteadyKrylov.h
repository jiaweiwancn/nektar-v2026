// SPDX-License-Identifier: MIT
#ifndef NEKTAR_WJ_NATIVE_STEADY_KRYLOV_H
#define NEKTAR_WJ_NATIVE_STEADY_KRYLOV_H

#include <Eigen/Core>
#include <array>
#include <functional>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace Nektar::WJNative
{
using SteadyVector=Eigen::VectorXd;

struct SteadyMeritState
{
    static constexpr const char *Identity="WJH27OriginalMeritV1";
    double usedMass=0., nextMass=0., referenceNorm=0., previousNorm=0., acceptedNorm=0.;
    double usedForcing=0., nextForcing=0., step=0.;
    std::optional<double> ratio;
    int accepted=0, goodSteps=0;
    std::string decision;
    std::string ToJson()const;
};

struct SteadyCurvatureOptions
{
    int maximumProducts=0;
    double maximumScaledSize=0., budgetFraction=0.;
};

struct SteadyMeritOptions
{
    double initialMass=0., reanchorMass=0.;
    int maximumWork=96;
    std::optional<SteadyMeritState> restored;
    double initialTrialFraction=0.;
    std::optional<SteadyCurvatureOptions> curvature;
};

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

// A bounded numerical preparation failure. The solver returns its last
// evaluated base and the ledger rather than treating it as a physical trial.
class SteadyCorrectionFailure : public std::runtime_error
{
public:
    using std::runtime_error::runtime_error;
};

struct SteadyTrialValidity
{
    bool admissible=true;
    std::string reason;
};

struct SteadyMeritRecord;
struct SteadyKrylovActionAudit;

struct SteadyKrylovCallbacks
{
    // All callbacks use PHYSICAL coordinates/residuals/corrections. Residual
    // must be the caller's current fixed equations with essential Gamma
    // eliminated. Steady preparation uses alpha=0; an implicit physical step
    // supplies its exact BDF mass/history. This engine adds no mass or time.
    std::function<SteadyVector(const SteadyVector &)> Residual;
    // Optional exact differentiation of a part of Residual. Remove its
    // current values from the same FD stencil, then add its base tangent.
    // The nonlinear residual, auxiliary loads and inverse remain separate.
    std::function<SteadyVector(const SteadyVector &)> ExactPartResidual;
    std::function<SteadyVector(const SteadyVector &,const SteadyVector &)> ExactPartAction;
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
    // Prepare the additive operator independently of the inverse, including
    // whenever a standalone directional probe supplies a new base.
    std::function<void(const SteadyVector &,const SteadyVector &)> PrepareCorrection;
    // Raw additive correction action in physical unknowns at the frozen
    // nonlinear base. It is used in every linear product and check, never F.
    std::function<SteadyVector(const SteadyVector &,const SteadyVector &)> CorrectionAction;
    // Used only by the opt-in merit policy. Set the mass at an already frozen
    // base before any linear product or inverse preparation.
    std::function<void(double)> SetCorrectionMass;
    std::function<void(const SteadyMeritState &)> ObserveMeritState;
    std::function<void(const SteadyMeritRecord &)> ObserveMeritDirection;
    // Opt-in only. The fixed unshifted colored action is a diagnostic control.
    std::function<void(const SteadyKrylovActionAudit&)> ObserveActionAudit;
    std::function<SteadyVector(const SteadyVector&)> ColoredAction, FactoredAction, TraceDirection;
    // For the merit policy TrueNorm is the maximum of these fixed quadratic
    // seminorms (including absolute value and identically zero blocks).
    std::function<std::array<double,6>(const SteadyVector &)> ResidualBlocks;
    // Optional infinity norms by free-variable block in scaled coordinates.
    std::function<std::array<double,6>(const SteadyVector &)> ScaledUnknownBlocks;
    // Diagnostic-only observations of the unshifted local model.
    std::function<void(int,const SteadyVector &,const SteadyVector &)> InspectDirection;
    std::function<void(int,int,double,const SteadyVector &,const SteadyVector &)> InspectTrial;
    std::function<void(const SteadyVector &,double)> ObserveAccepted;
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
    bool evaluateInitialOnly=false;
    int maximumInverseAge=0, reuseLinearIterations=12;
    bool measureDirection=false;
    bool adaptiveTrialStep=false;
    int frozenSizeColumns=0;
    std::optional<SteadyMeritOptions> merit;
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
    double plusDisplacement=0, minusDisplacement=0, halfDisplacement=0;
    double scaledDirectionInfinityNorm=0, seconds=0;
    // Admissible means ValidTrial + Residual domain checks only. These flags
    // do not certify optional full nonlinear-state/principal validation.
    bool plusAdmissible=false, minusAdmissible=false, halfAdmissible=false;
};

// Diagnostic payloads are separate from states and are never accepted iterates.
// Matrices use Eigen column-major order. The publisher binds the free ordering.
struct SteadyKrylovActionAudit
{
    std::string stage;
    int newton=-1, cycle=-1, iteration=-1, proposal=-1;
    long long factorSerial=0;
    std::vector<std::pair<std::string,double>> measurements;
    std::vector<std::pair<std::string,std::string>> notes;
    std::vector<std::pair<std::string,Eigen::MatrixXd>> arrays;
    std::vector<SteadyFDRecord> stencils;
    std::string FileSuffix(std::size_t sequence) const;
};

struct SteadyColoredColumnProbe
{
    int column=0, colorColumns=0;
    SteadyVector colorDirection, exactColorAction, exactSingleAction;
    std::vector<int> support;
};

struct SteadyLinearRecord
{
    int newton=0, cycle=0, iteration=0, leastSquaresRank=0;
    double projectedNorm=0, trueScaledNorm=0, trueNorm=0;
    std::optional<std::array<double,6>> trueBlocks;
    bool trueChecked=false, breakdown=false, converged=false, auxiliary=false;
};

struct SteadyLineSearchRecord
{
    int newton=0, trial=0, proposal=-1;
    double step=0, trueNorm=0;
    bool admissible=false, evaluated=false, sufficientDecrease=false, accepted=false;
    bool nonlinearValidationAttempted=false, nonlinearValidated=false;
    bool curvatureCorrected=false;
    std::optional<double> predictedNorm,decreaseRatio;
    std::optional<std::array<double,6>> predictedBlocks,actualBlocks,remainderBlocks;
    double scaledDisplacement=0.;
    std::optional<std::array<double,6>> displacementBlocks;
    std::string reason;
};

struct SteadyInverseUse
{
    int base=0,factorBase=0,age=0,iterations=0;
    bool rebuild=false,prepared=false;
    std::optional<bool> converged;
    std::optional<double> trueNorm;
    double preparationSeconds=0., linearSeconds=0.;
    long long serial=0;
    std::optional<double> mass;
    std::string reason;
};

struct SteadyMeritRecord
{
    int base=0, proposal=0, inverse=-1, work=0;
    double mass=0., forcing=0., baseNorm=0.;
    double linearNorm=0., scaledError=0., scaledBaseNorm=0.;
    double actionUncertainty=0., scaledUncertainty=0.;
    double scaledDirectionNorm=0., primaryParameter=0., halfParameter=0.;
    std::optional<std::array<double,6>> directionBlocks;
    std::optional<std::array<double,6>> uncertaintyBlocks;
    double step=0., predictedNorm=0., activeSlope=0., slopeUncertainty=0.;
    std::array<double,6> baseBlocks{}, shiftedErrorBlocks{}, correctionBlocks{},
        fullPredictedBlocks{}, predictedBlocks{};
    std::vector<int> activeBlocks;
    std::string decision,massMove;
    bool modelEvaluated=false, linearPassed=false;
};

struct SteadyCurvatureRecord
{
    int base=0, proposal=0, trial=0, products=0;
    long long factorSerial=0;
    double step=0., plainNorm=0., armijoLimit=0., primaryPrediction=0., seconds=0.;
    double scaledSizeLimit=0., scaledBudget=0.;
    std::array<double,6> defectBlocks{}, primaryUncertainty{}, blockBudgets{};
    std::optional<double> scaledCorrection, scaledError, scaledUncertainty,
        correctedNorm, correctedPrediction, correctedDisplacement;
    std::optional<std::array<double,6>> errorBlocks, uncertaintyBlocks, correctedBlocks,
        remainingDefectBlocks;
    std::string decision, reason;
};

struct SteadyNewtonRecord
{
    int iteration=0;
    double baseNorm=0, forcing=0, linearTrueNorm=0, acceptedNorm=0, acceptedStep=0;
    double acceptedSeconds=0.;
    bool linearEvaluated=false,linearConverged=false, accepted=false;
    std::string missingLinearReason;
};

struct SteadyKrylovResult
{
    bool converged=false, hasEvaluatedState=false, bestFromAcceptedIterate=true;
    bool nonlinearValidationEnabled=false;
    int bestNewtonIteration=-1,lastAcceptedNewtonIteration=-1;
    std::string status="NOT_CONVERGED", reason;
    // Best evaluated admissible state, with optional full nonlinear validation
    // passed, not necessarily an Armijo-accepted
    // iterate on failure. The caller must explicitly install any returned x.
    // If the deadline expires before an initial norm can be evaluated,
    // hasEvaluatedState=false and trueNorm/residual are not evidence.
    SteadyVector x, residual, lastAcceptedX;
    double trueNorm=0, lastAcceptedNorm=0;
    SteadyKrylovStatistics statistics;
    std::vector<SteadyFDRecord> finiteDifferences;
    std::vector<SteadyLinearRecord> linear;
    std::vector<SteadyLineSearchRecord> lineSearch;
    std::vector<SteadyNewtonRecord> newton;
    std::vector<SteadyInverseUse> inverses;
    std::vector<SteadyMeritRecord> meritDirections;
    std::vector<SteadyCurvatureRecord> curvatureTrials;
    std::vector<SteadyMeritState> meritAcceptances;
    std::optional<SteadyMeritState> meritStart, meritEnd;
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
    void AuditDirection(const SteadyVector &base,const SteadyVector &physicalDirection)const;
    void AuditDirections(const SteadyVector &base,const Eigen::MatrixXd &physicalDirections)const;
    void AuditColumn(const SteadyVector &base,const SteadyVector &physicalDirection,
                     const SteadyColoredColumnProbe &column,const std::string &label)const;
private:
    SteadyKrylovCallbacks m_callbacks;
    SteadyVector m_unknownScale,m_residualRowScale;
    SteadyKrylovOptions m_options;
};
}
#endif
