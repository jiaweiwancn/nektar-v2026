// SPDX-License-Identifier: MIT
#ifndef NEKTAR_WJ_NATIVE_COUPLING_H
#define NEKTAR_WJ_NATIVE_COUPLING_H
#include <IncNavierStokesSolver/EquationSystems/WJNativeMixedFlow.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeScalar.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeReference.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeCoupledMMS.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeWallV2.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeSteadyKrylov.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeSteadyCorrection.h>
#include <filesystem>
#include <algorithm>
#include <map>
#include <optional>

namespace Nektar::WJNative
{
struct CoupledState { Vec flow,k,epsilon; };
struct CoupledHistory
{
    std::uint64_t step=0;
    double time=0,dt=0;
    CoupledState current,previous,advection,previousAdvection;
};
struct CoupledRecoveryTiming
{
    double readSeconds = 0;
    double cacheRebuildSeconds = 0;
};
struct CoupledOptions
{
    double dt=.005,linearTolerance=1e-12,nonlinearTolerance=1e-10,mmsBodyAmplitude=0;
    int maximumIterations=60,quadratureX=8,quadratureR=12,andersonDepth=0;
    bool manufactured=false,recordNonlinearTrace=false,checkCurrentActions=false;
    double hillInitialSlowdown=0., hillInitialEpsilonBump=0., geometrySeedTaper=1.;
    double geometrySeedFootLength=0., geometrySeedCollarLength=0., geometrySeedCollarReach=0.;
    std::string geometrySeedPacketSha256;
    std::string outletPolicy="FullTractionV1";
    bool checkOutletActions=false;
    double diagnosticOutletX=15.;
    // Numerical split only; the complete current V is identical in both.
    std::string stressSplit="PositiveShearV1",steadyVelocityPC="FullGradientV1";
    std::string advancement="PicardV1";
    std::string flowPreconditioner="SparseLU", picardStop="CoupledResidual";
    std::string picardAdvection="Explicit", picardStressLinearization="FrozenStrain";
    std::string newtonAdvection="Explicit";
    int pressureCorrectors=1, picardOuterCorrectors=0;
    SteadyKrylovOptions stepNewton;
    std::string stepPreconditioner="RebuildEveryNewtonV1";
    double stepPCRelativeDrift=.02;
    int stepPCMaximumUses=8;
    ConstitutiveBackend backend=ConstitutiveBackend::Historical;
    std::string sourceIdentity,mmsVariant="Historical",initialSourceIdentity,initialReferenceIdentity;

    bool UsesImplicitAdvection() const
    {
        return advancement == "CoupledJFNKV1" ? newtonAdvection == "Implicit"
                                              : picardAdvection == "Implicit";
    }
};
struct CoupledResidual
{
    double momentum=0,continuity=0,gauge=0,k=0,epsilon=0,trace=0;
    double Maximum()const{return std::max({momentum,continuity,gauge,k,epsilon,trace});}
};
enum class DiagnosticRole { RejectedInitialization, RejectedTimeStep, SteadySearch };
struct DiagnosticCandidate
{
    DiagnosticRole role;
    CoupledState state;
    std::uint64_t attemptedStep=0;
    double attemptedTime=0.;
    std::string reason;
    std::optional<CoupledResidual> residual;
    int attemptedIteration=-1;
};
struct DiagnosticSteadyOptions
{
    // Zero pseudoTime disables the correction mass. These shifts never enter F.
    double pseudoTime=0., correctionDiffusion=0.;
    // A positive eta selects a separate momentum problem for initialization.
    double eta=0., regularizationViscosity=0.;
    double regularizationStart=0., regularizationFull=0.;
    std::string seedSha256;
    int scalarSubdivisionDepth=0;
    bool sparseInverse=false, auditCorrection=false;
};
struct CoupledTiming
{
    double constitutive=0,wall=0,assembly=0,flow=0,scalar=0,rebuild=0,
           guards=0,forcing=0,residual=0,nonlinear=0,io=0,other=0;
    double Sum()const{return constitutive+wall+assembly+flow+scalar+rebuild+guards+forcing+residual+nonlinear+io+other;}
};
struct CoupledStepTimingDetail
{
    double momentumAdvection = 0;
    double scalarAdvection = 0;
    double scalarMomentumCoupling = 0;
    double historyMovement = 0;
    double kFactor = 0;
    double epsilonFactor = 0;
};
struct CoupledStepInfo
{
    int nonlinearIterations=0,flowIterations=0,kIterations=0,epsilonIterations=0,rebuilds=0,retries=0;
    std::string retryReason;
    double seconds=0,coefficientDrift=0,minimumK=0,minimumEpsilon=0,
           minimumCovarianceOverK=0,minimumPrincipalOverMu=0;
    CoupledResidual residual;
    CoupledTiming timing;
    CoupledStepTimingDetail timingDetail;
    MixedFlowWork flowWork;
    bool acceptedByFixedBudget = false;
    WallEvaluationStats constitutiveStats;
    struct Iteration {CoupledResidual residual;bool accelerated=false;double seconds=0;};
    std::vector<Iteration>nonlinearTrace;
    int accelerationProposals=0,accelerationAccepted=0,accelerationInadmissible=0;
    long long coupledKrylovIterations=0,coupledJacobianActions=0,
              coupledResidualEvaluations=0,coupledPreconditionerApplications=0,
              coupledRejectedTrials=0;
    bool coupledStatisticsComplete=true;
    int timePCRefreshes=0,timePCReuses=0,timePCMaximumUsedAge=0;
    // Exclusive reasons: cold, alpha, drift, age, control, retry, injected.
    std::array<int,7>timePCRefreshReasons{};
    struct TimePCEvent
    {
        int attempt=0,base=0,ageBefore=0,ageUsed=0;
        double alpha=0,drift=0;
        std::string reason;
        bool refreshed=false;
    };
    std::vector<TimePCEvent>timePCEvents;
};

// Append comma-prefixed numerical cache members to an open diagnostic object.
void WriteStepPreconditionerFields(std::ostream&out,const CoupledOptions&options,
                                  const CoupledStepInfo&info);
void WriteFlowSolverFields(std::ostream &out, const CoupledOptions &options,
                           const CoupledStepInfo &info);

/** Serial coupled BE/BDF2 velocity/modified-pressure,
 * full-domain k and outer-only epsilon. E_w and the reference cutoff are
 * current dependent data. No artificial equilibrium source or scalar floor.
 */
class Coupling
{
public:
    Coupling(std::vector<MultiRegions::ContFieldSharedPtr> fields,Geometry geometry,
             std::shared_ptr<const Reference> reference,CoupledOptions options);
    CoupledHistory Initialise()const;
    CoupledHistory FreshHistory(const CoupledState&state)const;
    CoupledHistory ReadInitialFields(const std::string &path, bool allowCandidate=false,
                                     std::optional<double> geometrySeedHeight=std::nullopt);
    CoupledState PrepareState(const std::string &path, const std::string &seedSha256,
                              std::optional<double> geometrySeedHeight);
    void EnableRejectionPackets(bool enabled){m_exportRejections=enabled;}
    void WriteRejection(const std::filesystem::path&path)const;
    SteadyKrylovResult SolveDiagnosticSteady(const CoupledState &seed,
        const SteadyKrylovOptions &options,CoupledStepInfo &final,double preconditionerMass,
        const DiagnosticSteadyOptions &diagnostic);
    void WriteSteadyDiagnostic(const std::filesystem::path &path)const;
    CoupledHistory PerturbInitialHistory(const CoupledHistory&base,const std::string&certificate,
                                        const std::string&sha256);
    CoupledHistory SteadyCandidateHistory(const Vec&unknowns)const;
    CoupledStepInfo AuditSteadyState(const CoupledState&state,double time=0);
    CoupledStepInfo Advance(CoupledHistory&history,bool injectRetry=false);
    const CoupledStepInfo&LastAttemptInfo()const{return m_lastAttempt;}
    CoupledResidual SteadyResidual(const CoupledState&state,double time=0);
    /** Native steady equations with dependent Gamma eliminated. A failed
     * reference solve preserves the caller's complete history. On success a
     * fresh step-zero history is constructed; Newton iterations are never
     * interpreted as elapsed physical time. */
    SteadyKrylovResult SolveSteady(CoupledHistory&history,const SteadyKrylovOptions&options,
                                 CoupledStepInfo&final,double preconditionerMass=0);
    const MixedFlow&FlowKernel()const{return m_flow;}
    const Scalar&KKernel()const{return m_k;}
    const Scalar&EpsilonKernel()const{return m_e;}
    const RuleFactory&Factory()const{return m_rules;}
    const Reference&ReferenceModel()const{return *m_reference;}
    const CoupledMMS&ManufacturedSolution()const{return *m_mms;}
    const CoupledOptions&Options()const{return m_options;}
    const std::string&Signature()const{return m_signature;}
    double TopRadius()const{return m_top;}
    void Store(const CoupledState&state)const;
    void WriteHistory(const std::filesystem::path&path,const CoupledHistory&history)const;
    CoupledHistory ReadHistory(const std::filesystem::path &path,
                               CoupledRecoveryTiming *recoveryTiming = nullptr) const;
    void Report(std::ostream&out,const CoupledHistory&history,const CoupledStepInfo&info,
                bool includeCoefficients=false)const;
    std::array<double,5> MMSErrors(const CoupledState&state,double time)const;
    void ResetPreconditioners();
private:
    friend struct InnerKFixtureAccess;
    CoupledStepInfo AdvanceNewton(CoupledHistory&history,bool injectRetry);
    CoupledState InitialState()const;
    CoupledState LiftInitialState(const CoupledState&state)const;
    CoupledState ImportInitialState(const std::string &path,bool allowCandidate,
        std::optional<double> geometrySeedHeight,const std::string *diagnosticSeedSha256);
    SteadyKrylovResult SearchSteady(const CoupledState &seed,
        const SteadyKrylovOptions &options,CoupledStepInfo &final,double preconditionerMass,
        CoupledHistory *physicalHistory,const DiagnosticSteadyOptions *diagnostic);
    struct Data
    {
        struct VelocityInput {double k,E,aux,r;Tensor gradient;Vector normal;};
        std::vector<std::vector<VelocityInput>>velocityInputs;
        Rules kRules;
        std::vector<Vec>mu,bulk,shear,Dk,De,sourceK,sourceE,reactionE;
        std::vector<std::vector<Tensor>>remainder,viscousStress;
        Vec coefficients,kCoefficients,eCoefficients;
        double minK=1e100,minE=1e100,minCovariance=1e100,minPrincipal=1e100,
               minModifiedPrincipal=1e100,traceError=0;
    };
    struct Loads { Vec flow,k,epsilon; };
    void BuildTopology();
    const Element&Owner(int geometry)const;
    const Element&HOwner(const Element&element)const;
    NativeRule HRule(const NativeRule&rule)const;
    const Element& MatchingOwner(double a) const;
    double InterfaceK(const Element&element,double x,const CoupledState&state)const;
    std::array<double,2> InterfaceBounds(const CoupledState &state,bool requireQualified=true) const;
    std::array<double,2> InterfaceVelocity(const Element&element,double x,const CoupledState&state)const;
    std::vector<double> CutoffXBreaks(const Element&element,const std::array<double,4>&coefficients,
                                    double traceLeft,double traceRight)const;
    std::vector<double> CurrentCutoffXBreaks(const Element&element,const CoupledState&state)const;
    std::vector<double> MMSCutoffXBreaks(const Element&element,double time)const;
    void CheckTraceEvaluation()const;
    double NativeRadius(double physicalZ)const;
    double NativeX(double x)const;
    CoupledState DirichletState(const CoupledState&state,double time)const;
    enum class PrincipalPolicy { Skip, Observe, Enforce };
    Data Evaluate(const CoupledState&state,double time,PrincipalPolicy principal,CoupledTiming&timing,
                  const DiagnosticSteadyOptions *diagnostic=nullptr)const;
    Data Evaluate(const CoupledState&state,double time,bool principal,CoupledTiming&timing)const
    {return Evaluate(state,time,principal?PrincipalPolicy::Enforce:PrincipalPolicy::Skip,timing);}
    Loads Forcing(double time,CoupledTiming&timing,bool freezeManufactured=false)const;
    double RegularizationShape(const NativeRule &rule,int q,const DiagnosticSteadyOptions &options)const;
    Sparse Strain(const Data&data)const;
    Sparse VelocityTangent(const Data&data,CoupledTiming&timing)const;
    Vec TimePreconditionerData(const Data&data)const;
    Vec RemainderLoad(const Data&data)const;
    Loads CurrentActions(const CoupledState&state,const Data&data)const;
    Vec MomentumAdvection(const CoupledState&state)const;
    Sparse MomentumTransportJacobian(const CoupledState&state)const;
    CoupledState Advection(const CoupledState &state, const Rules &kRules,
                           CoupledStepTimingDetail *timing = nullptr) const;
    CoupledResidual Residual(const CoupledState&state,const Data&data,const Loads&base,
                             double alpha,double time,CoupledTiming&timing,
                             const CoupledHistory*history=nullptr)const;
    Vec FlowLift(double time)const;
    Vec DependentTraceDirection(const CoupledState &state, const Vec &dk) const;
    Vector InitialVelocity(double a, double r) const;
    struct SeedGradientFit
    {
        Sparse gram;
        Vec load;
        double targetSquared=0., weight=0.;
        double Error(const Vec &velocity) const;
    };
    SeedGradientFit GeometrySeedGradientFit(const Vec &oldFlow, const Geometry &source) const;
    Vec ProjectContinuity(const Vec &flow, const SeedGradientFit *fit=nullptr) const;
    double ContinuityNorm(const Vec &flow) const;
    Vec RemapVelocity(const Vec &flow, const Geometry &source);
    void WriteOutletActions(std::ostream&out,const CoupledHistory&history)const;
    Vec PackFree(const CoupledState&state)const;
    CoupledState FromFree(const Vec&free,double time)const;
    Vec UnknownScales()const;
    std::unique_ptr<SteadySparseInverse> MakeSteadySparseInverse(const Vec &unknownScale,const Vec &rowScale)const;
    Vec KVolumeLoad(const CoupledState&state,const Data&data)const;
    Vec EVolumeLoad(const CoupledState&state,const Data&data)const;
    double Drift(const Vec&current,const Vec&previous)const;
    void CheckHistory(const CoupledHistory&history)const;
    void WriteCoefficientFields(std::ostream&out,const CoupledState&state,double time,bool refreshLift=true)const;
    void WriteInputHistory(std::ostream&out)const;
    void WriteDiagnostic(const std::filesystem::path&path,const DiagnosticCandidate&candidate,
                         const CoupledStepInfo *original,const DiagnosticSteadyOptions *options)const;
    std::vector<MultiRegions::ContFieldSharedPtr>m_fields;
    Geometry m_geometry;
    std::shared_ptr<const Reference>m_reference;
    CoupledOptions m_options;
    MixedFlow m_flow;
    RuleFactory m_rules;
    Scalar m_k,m_e;
    std::unique_ptr<CoupledMMS>m_mms;
    double m_left=0,m_right=0,m_top=0,m_volume=0;
    std::map<int,int>m_geometryIndex,m_hOwner;
    struct KTrace {std::array<int,4>ids;std::array<double,4>signs;};
    std::map<int,KTrace>m_hKTrace;
    std::array<double,2>m_cutoffRange{};
    std::map<double,double>m_cutoffFaceK;
    std::vector<std::pair<int,int>>m_boundaryFaces;
    std::string m_signature,m_outletIdentity;
    bool m_zeroTransverseOutlet=false;
    std::map<int,bool>m_outletReactionRows; // Vertical coefficient ID and shared-corner flag.
    Vec m_preconditionerCoefficients;
    Vec m_kPreconditionerCoefficients,m_ePreconditionerCoefficients;
    double m_flowFactorAlpha=-1,m_kFactorAlpha=-1,m_eFactorAlpha=-1;
    // These are inverse approximations only. Current residual/Jv/constitutive
    // data never read them, and no derived factor is serialized in a restart.
    std::array<Sparse,3>m_timePCOperators;
    Sparse m_timePCEpsilonFull;
    Vec m_timePCData;
    double m_timePCAlpha=-1;
    int m_timePCUses=0;
    bool m_timePCReady=false;
    CoupledHistory m_lastInputHistory;
    CoupledStepInfo m_lastAttempt;
    bool m_hasInputHistory=false;
    bool m_exportRejections=false;
    mutable std::optional<DiagnosticCandidate> m_rejectedCandidate;
    std::optional<DiagnosticCandidate> m_steadyDiagnosticCandidate;
    DiagnosticSteadyOptions m_diagnosticOptions;
    CoupledStepInfo m_diagnosticOriginal;
    double m_diagnosticExtraWork=0.;
    double m_diagnosticModifiedMargin=0.;
    std::array<double,2> m_diagnosticWholeCellLower{};
    std::vector<std::array<double,3>> m_diagnosticCorrections;
    std::vector<SteadySparseWork> m_sparseSteadyWork;
    std::vector<std::array<double,6>> m_steadyBaseObservations;
    struct SteadyActionAudit
    {
        std::string direction, rejection;
        double step=0., taylorRemainder=0., centeredRelativeError=0., shiftRelativeError=0.;
        std::optional<std::array<double,6>> differenceBlocks;
    };
    std::vector<SteadyActionAudit> m_steadyActionAudits;
    std::array<std::array<double,4>,4> m_steadyScaleExtrema{};
    double m_steadyPatternSeconds=0., m_steadyAuditSeconds=0.;
    int m_steadyAuditResiduals=0;
    std::string m_steadyBaseSha256;
    mutable Vec m_bulkBodyAction,m_sampledBodyAction;
    mutable Vec m_outletLoad;
    mutable Vec m_physicalFlowLift,m_physicalKLift;
    mutable bool m_hasPhysicalForcing=false;
    mutable Loads m_physicalForcing;
    std::string m_initialPacketSha,m_initialSourceStateSha;
    bool m_initialProjected=false,m_initialCandidate=false;
    std::optional<double> m_geometrySeedHeight;
    bool m_geometrySeedRefined=false;
    std::uint64_t m_geometrySeedStep=0;
    double m_geometrySeedTime=0.;
    double m_initialFluxRoundtripError=0.,m_initialTangentialTraceCorrection=0.;
    std::optional<std::array<double,2>> m_initialGradientFit;
    mutable WallEvaluationStats m_constitutiveStats;
    mutable std::uint64_t m_actionChecks=0,m_augmentedChecks=0;
    mutable std::array<double,3>m_actionDifference{};
    mutable double m_augmentedDifference=0,m_checkedWeakDivergence=0,m_checkedGaugeMultiplier=0;
    mutable double m_dependentTraceDerivativeError=0.;
    mutable double m_initialContinuityBefore=0.,m_initialContinuityAfter=0.;
};
}
#endif
