///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_RANS_WALL_EVOLUTION_PATCH_H
#define NEKTAR_RANS_WALL_EVOLUTION_PATCH_H
#include <IncNavierStokesSolver/EquationSystems/RANSWallGlobalField.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSharedInventory.h>
#include <LibUtilities/Communication/Comm.h>
#include <SolverUtils/Turbulence/RANSWallNormalMap.h>
#include <SolverUtils/Turbulence/EquilibriumAWFTransition.h>
#include <SolverUtils/Turbulence/RANSRobustVariant.h>
#include <SolverUtils/Turbulence/RANSWallInnerStress.h>
#include <algorithm>
#include <cstdint>
#include <set>

namespace Nektar
{
namespace RANSWallR3 { class DurableProgress; }
SolverUtils::RANSWallNormalMap LoadRANSWallNormalMap(const LibUtilities::SessionReaderSharedPtr &session,int dimension,double height);
using RANSWallEvolutionFunction=std::function<double(const RANSWallCoordinate&,double)>;
struct RANSWallEvolutionData
{
    std::vector<RANSWallEvolutionFunction> innerBody,outerBody,boundary;
    // Physical exterior mean tractions and outward diffusive scalar fluxes.
    // Session functions/data are immutable within generation; Touch on changes.
    std::vector<RANSWallEvolutionFunction> outletInner,outletOuter,topTraction;
    RANSWallEvolutionFunction topKFlux,topEpsilonFlux;
    // Sensors retain the physical body force when an IMEX stage moves its
    // prescribed load to the explicit side. Empty means use outerBody.
    std::vector<RANSWallEvolutionFunction> sensorBody;
    RANSWallEvolutionFunction innerKSource,outerKSource,outerEpsilonSource;
    // Optional fixed native external traces, e.g. the stored neutral target.
    std::map<int,double> fixedBoundary;
    // Callers changing a callable or its captured data inside an evaluation
    // scope must advance this immutable-data generation. Scope-local reuse
    // also compares fixed-boundary bytes and exact physical time.
    std::uint64_t generation=0;
    void Touch(){++generation;}
    explicit RANSWallEvolutionData(int dimension);
};

/** Analytically mapped 2D/3D native CG stage with globally ordered replicated coefficients.
 * Physical cell quadrature is owned once across the supplied communicator;
 * residuals are collectively assembled and one root solves the dense stage.
 * This bounded reference backend is deliberately distinct from a scalable
 * distributed sparse solver. Inner/outer pressure and all histories survive.
 */
class RANSWallEvolutionPatch
{
    struct EvaluationCache;
public:
    using Vector=RANSWallAlgebra::Vector;
    using WallMatrix=RANSWallAlgebra::DenseMatrix;
    using BlockMatrix=RANSWallAlgebra::MappedMatrix;
    using Diagnostic=std::map<std::string,Vector>;
    class EvaluationScope
    {
        const RANSWallEvolutionPatch &patch;
        std::shared_ptr<EvaluationCache> previous;
        bool owns=false;
    public:
        EvaluationScope(const RANSWallEvolutionPatch &,const RANSWallEvolutionData &,double);
        ~EvaluationScope();
        EvaluationScope(const EvaluationScope&)=delete;
    };
    std::string EvaluationReuseJSON() const;
    std::string EvaluationReuseChecks(const Vector &,const Vector &,double,const RANSWallEvolutionData &);
    RANSWallEvolutionPatch(const LibUtilities::SessionReaderSharedPtr &serialSession,
        const SpatialDomains::MeshGraphSharedPtr &graph,
        const LibUtilities::CommSharedPtr &collective,const RANSWallPatchParameters &parameters,
        const SolverUtils::RANSWallNormalMap &mapping=SolverUtils::RANSWallNormalMap{});
    int Dimension() const { return m_dimension; }
    int K() const { return m_dimension; }
    int Epsilon() const { return m_dimension+1; }
    int InnerPressure() const { return m_dimension+2; }
    int OuterPressure() const { return m_dimension+3; }
    int Gauge() const { return m_dimension+4; }
    int DifferentialSize() const { return offset[InnerPressure()]; }
    int Size() const { return offset[Gauge()]+size[Gauge()]; }
    bool HasGauge() const { return size[Gauge()]!=0; }
    double GaugeValue(const Vector &s) const { return HasGauge()?s[offset[Gauge()]]:0.; }
    bool OpenBoundary() const { return m_boundaryMode=="OpenTunnelV1"; }
    const std::string &InitialVelocityMode() const { return m_initialVelocityMode; }
    int InitialVelocityCalls() const { return m_initialVelocityCalls; }
    static void ValidateFlatInletTraceMap(const SolverUtils::RANSWallNormalMap &);
    Vector ApplyInitialVelocity(const Vector &state);
    void LoadOpenBoundaryData(const LibUtilities::SessionReaderSharedPtr &,RANSWallEvolutionData &) const;
    int PhysicalAxis(int a) const { return m_dimension==2 && a==1?2:a; }
    const RANSWallPatchParameters &Parameters() const { return m_parameters; }
    const LibUtilities::CommSharedPtr &Comm() const { return m_comm; }
    const SolverUtils::RANSWallNormalMap &Mapping() const { return m_mapping; }
    bool Mapped() const { return !m_mapping.Identity(); }
    bool MappedBlock() const { return m_linearBackend=="MappedBlockV1"; }
    int Owner(int geometry) const { return (geometry+m_ownershipShift)%m_comm->GetSize(); }
    RANSWallCoordinate Physical(const RANSWallCoordinate &chart) const;
    RANSWallGlobalPoint Point(int block,const RANSWallCoordinate &chart,int side=0) const;
    std::vector<RANSWallGlobalPoint> Quadrature(int block,int order=0) const;
    struct FacePoint { RANSWallGlobalPoint test;SolverUtils::RANSVector area;double weight; };
    struct ExteriorPoint { FacePoint face;std::string role;bool inner; };
    std::vector<ExteriorPoint> Exterior(int order=0) const;
    Vector ExteriorResidual(const Vector &,double,const RANSWallEvolutionData &,int region=-1) const;
    Vector ExteriorPressureAction(const Vector &,int region=-1) const;
    void ValidateOutlet(const Vector &,const std::string &phase) const;
    std::map<int,double> OpenConstraints(double,const RANSWallEvolutionData &) const;
    std::string OpenBoundaryDiagnostic(const Vector &,const Vector &,double,double,const RANSWallEvolutionData &,
        const std::vector<RANSWallCoordinate> &physicalLocations={}) const;
    std::vector<FacePoint> Surface(int block,int axis,double position,int side,int order=0) const;
    double Distance(const RANSWallCoordinate &chart) const;
    double DependentEpsilon(const Vector &state,const RANSWallCoordinate &chart) const;
    double DependentEpsilon(const Vector &state,const RANSWallCoordinate &chart,double time,
        const RANSWallEvolutionData &data) const;
    bool UniformEPFlat() const { return m_uniformEP; }
    bool WJ() const { return m_parameters.model.Model()==SolverUtils::RANSModel::WJHighRe; }
    bool WJTerrain() const { return WJ() && m_nativeGeometry=="BudgetTerrainV1"; }
    const SolverUtils::WJParameters &WJConfiguration() const { return m_parameters.wj; }
    std::string EpsilonName() const { return WJ()?"epsilon":"eps_tilde"; }
    double SigmaK() const { return WJ()?m_parameters.wj.sigmaK:m_parameters.model.Parameters().sigmaK; }
    double Epsilon1() const { return WJ()?m_parameters.wj.cEpsilon1:m_parameters.model.Parameters().cEpsilon1; }
    double Epsilon2() const { return WJ()?m_parameters.wj.cEpsilon2:m_parameters.model.Parameters().cEpsilon2; }
    bool BudgetAWF() const { return m_formulation=="CLSCalibratedBudgetAWF" || m_formulation=="WJCalibratedBudgetAWF"; }
    bool AcceptedAWF() const { return BudgetAWF(); } // Legacy compatibility.
    bool WJEP() const { return WJ() && m_formulation=="WJCalibratedEPAWF"; }
    bool UncappedWall() const { return !m_cappedDestruction; }
    const std::string &FormulationName() const { return m_formulation; }
    bool MatchingEnabled() const { return m_transitionEnabled || BudgetAWF() || WJEP(); }
    std::string WallPolicyName() const { return (BudgetAWF() || WJEP())?m_formulation:(m_uniformEP?"UniformEPFlat-v1":(m_transitionEnabled?"Upstream":"Fixed")); }
    bool TransitionEnabled() const { return m_transitionEnabled; }
    const SolverUtils::AWFTransitionControls &TransitionControls() const { return m_transition; }
    const Vector &TransitionReference() const { return m_transitionReference; }
    std::array<double,3> TransitionDistances() const { return {m_sensorA,m_sensorB,m_baselineStart}; }
    double WallCoordinate(const RANSWallCoordinate &chart) const;
    SolverUtils::AWFTransitionInput TransitionInput(const Vector &state,const RANSWallCoordinate &chart,
        double time,const RANSWallEvolutionData &data) const;
    SolverUtils::AWFTransitionValues Transition(const Vector &state,const RANSWallCoordinate &chart,
        double time,const RANSWallEvolutionData &data) const;
    std::shared_ptr<SolverUtils::EquilibriumAWFProfiles> StateProfile(const Vector &state,
        const RANSWallCoordinate &chart,double time,const RANSWallEvolutionData &data) const;
    // Baseline and correction are projected separately. The second map is
    // constrained to zero on every designated baseline face and shared mode.
    std::array<std::map<int,double>,2> DissipationTrace(const Vector &state,double time,
        const RANSWallEvolutionData &data) const;
    std::array<SolverUtils::RANSVector,2> ForcingParts(const Vector &state,
        const RANSWallCoordinate &chart,double time,const RANSWallEvolutionData &data) const;
    /** Exact discrete affine pressure/gauge action. Differential entries of
     * direction must be zero. Supplied strong rows are masked, as in Residual.
     * This avoids subtracting pressure-independent stage loads. */
    Vector AffinePressureAction(const Vector &direction,
        const std::map<int,double> &constraints={},int region=-1) const;
    /** Complete pressure-only directional derivative, including transition
     * sensors, capped source normalization and projected dependent trace.
     * History, alpha and the fixed anchor have zero pressure derivative.
     * At a cap max-knot the existing directional kernel selects its actual
     * one-sided action; no two-sided smoothness is presumed there. */
    Vector PressureDirectionalAction(const Vector &state,const Vector &direction,
        double time,const RANSWallEvolutionData &data,bool imposeConstraints=true,
        int region=-1) const;
    Vector DirectionalAction(const Vector &state,const Vector &direction,const Vector &history,
        double alpha,double time,const RANSWallEvolutionData &data,double step=2e-6) const;
    /** Opt-in R3 spatial nonpressure actions, including moving normalization
     * rules and projected dependent-trace constraints. Time mass and the
     * separately qualified pressure/gauge columns are supplied by the caller.
     */
    WallMatrix NonpressureJacobianR3(const Vector &state,double time,
        const RANSWallEvolutionData &data,std::vector<std::string> *metadata=nullptr) const;
    Vector NonpressureDirectionalActionR3(const Vector &state,const Vector &direction,
        double time,const RANSWallEvolutionData &data,bool imposeConstraints=true,int region=-1) const;
    double BoundaryFluxCorrection() const { return m_boundaryCorrection; }
    std::vector<std::shared_ptr<RANSWallGlobalField>> fields;
    std::shared_ptr<RANSWallGlobalField> innerVelocity,outerVelocity;
    std::vector<int> offset,size;
    RANSWallCoordinate lower{},upper{};
    Vector Project(const std::vector<RANSWallEvolutionFunction> &functions,double time) const;
    Vector NeutralReference(double ustar,const std::array<double,2> &direction) const;
    std::map<int,double> FixedReferenceBoundary(const Vector &reference) const;
    std::map<int,double> Constraints(const Vector &state,double time,const RANSWallEvolutionData &data) const;
    Vector Lift(const Vector &state,double time,const RANSWallEvolutionData &data) const;
    Vector Residual(const Vector &state,const Vector &history,double alpha,double time,
        const RANSWallEvolutionData &data,bool anchored=true,bool imposeConstraints=true,int region=-1) const;
    Vector Solve(Vector state,const Vector &history,double alpha,double time,const RANSWallEvolutionData &data,
        bool condensed,int &iterations,double tolerance=2e-12,int maximumIterations=20) const;
    struct ScalarObservation
    {
        std::string inventory;
        RANSWallCoordinate coordinate{};
        int block=0,geometry=-1,side=0;
        double value=0,scale=1,weight=0;
        bool strict=false,dependentTrace=false;
        std::shared_ptr<RANSWallGlobalPoint> functional;
    };
    static const char *AuditRevision() { return "NativeElementSideV2"; }
    struct ObservationProvenance
    {
        std::shared_ptr<const RANSWallGlobalPoint> velocity;
        bool inner=false;
        int radialSide=0;
    };
    ObservationProvenance ResolveObservation(const ScalarObservation &point) const;
    RANSWallGlobalPoint GeometryPoint(int block,int geometry,const RANSWallCoordinate &coordinate) const;
    std::string InterfaceProvenanceChecks() const;
    struct ScalarDomainResult
    {
        bool admissible=false;
        bool normalizationApplicable=true;
        std::string phase,failure,limitingInventory,profileBranches;
        RANSWallCoordinate limitingPoint{};
        std::shared_ptr<RANSWallGlobalPoint> limitingFunctional;
        double minimumMargin=1e300,minimumInnerK=1e300,minimumOuterK=1e300,
            minimumOuterEpsilon=1e300,minimumSampledK=1e300,minimumTraceEpsilon=1e300,
            minimumB=1,minimumLambda=1,minimumCutoff=1e300,maximumCutoff=-1e300;
        double minimumBudgetRatio=1e300,maximumBudgetRelativeDefect=0.;
        std::map<std::string,int> counts;
        Vector effectiveState,traceState;
        RANSWallSharedInventory<ScalarObservation> observations;
        struct PendingProfile
        {
            RANSWallCoordinate anchor;
            std::shared_ptr<SolverUtils::EquilibriumAWFProfiles> profile;
            bool source=false;
        };
        std::vector<PendingProfile> pendingProfiles;
    };
    // Scalar evaluation domain only: covariance and the restricted outer
    // physical box are separate accepted-stage gates.
    ScalarDomainResult ScalarDomain(const Vector &state,double time,
        const RANSWallEvolutionData &data,bool retainLocations=false,bool kInventoryOnly=false,
        const std::vector<RANSWallCoordinate> &additionalAnchors={}) const;
    static std::string ScalarDomainJSON(const ScalarDomainResult &domain,bool locations=false);
    Diagnostic CompleteRestrictedAudit(const Vector &state,const Vector &history,
        double alpha,double time,const RANSWallEvolutionData &data,bool rangesOnly=false) const;
    WallMatrix ScalarFeasibleJacobian(const Vector &state,const Vector &history,
        double alpha,double time,const RANSWallEvolutionData &data,
        std::vector<std::string> *columnMetadata=nullptr,bool scalarOnly=false) const;
    Vector SolveFrozenScalars(Vector state,const Vector &history,double alpha,double time,
        const RANSWallEvolutionData &data,int &iterations,double tolerance,int maximumIterations) const;
    Vector SolveScalarFeasibleR3(Vector state,const Vector &history,double alpha,double time,
        const RANSWallEvolutionData &data,int &iterations,double tolerance,
        int maximumIterations,bool benchmarkOnly=false,bool localR4=false,bool localR5=false,
        const WallMatrix &frozenMatrix={},const std::string &frozenAssembly="[]") const;
    const std::string &R3ProgressFile() const { return m_r3ProgressFile; }
    const Vector &LastFeasibleState() const { return m_lastFeasibleState; }
    const std::string &NonlinearPolicyName() const { return m_nonlinearPolicy; }
    const std::string &ActionPolicyName() const { return m_actionPolicy; }
    bool CoupledR4() const { return m_globalization=="CoupledR4"; }
    bool JointRestorationR5() const { return m_globalization=="JointRestorationR5"; }
    bool DirectionReplayR6() const { return m_directionAccuracy=="ReplayR6"; }
    bool RefinedDirectionR6() const { return m_directionAccuracy=="RefinedR6"; }
    bool CompareR7() const { return m_r7Policy=="Compare"; }
    bool InteriorProbeR7() const { return m_r7Policy=="InteriorProbe"; }
    bool AdaptiveR7() const { return m_r7Policy=="Adaptive"; }
    bool GroundSafeguard() const { return CoupledR4() || JointRestorationR5(); }
    bool PimpleSingleSweep() const { return m_temporalCoupling=="PimpleSingleSweep"; }
    // Execution-only dense capacity. No constructor size restriction or
    // physical checkpoint identity change accompanies this finite opt-in.
    int DenseEvolutionCapacity(bool condensed) const;
    void CheckDenseEvolutionCapacity(bool condensed) const;
    std::string DenseEvolutionCapacityJSON() const;
    std::string DenseEvolutionPreflight(const Vector &,double,
        const RANSWallEvolutionData &,bool condensed,bool enforce) const;
    const std::string &TemporalCouplingName() const { return m_temporalCoupling; }
    Vector InitializationWork() const { return {double(m_liftCalls),double(m_liftIterations)}; }
    const std::vector<std::string> &PimpleRecords() const { return m_pimpleRecords; }
    Diagnostic PhysicalSamples(const Vector &state,double time,const RANSWallEvolutionData &data) const;
    bool PimpleKBoundEnabled() const { return m_pimpleKBound!="Off"; }
    static Vector ContinuationProfileAbscissae(){return {-3.,-1.,0.,1.,1.6};}
    Vector ComparisonProfileAbscissae() const
    {return ProspectiveKInventory()?Vector{-3.5,-2.625,-1.,1.6}:Vector{-3.5,-1.,1.6};}
    Vector ConstitutiveReportingAbscissae() const
    {
        auto values=ContinuationProfileAbscissae();
        if(ProspectiveKInventory())for(double x:ComparisonProfileAbscissae())
            if(std::find(values.begin(),values.end(),x)==values.end())values.push_back(x);
        return values;
    }
    bool ConstitutiveKInventory() const { return PimpleKBoundEnabled() &&
        (m_pimpleKInventory=="ConstitutiveV2" || m_pimpleKInventory=="ConstitutiveV3" || m_pimpleKInventory=="AcceptedBudgetV1" || m_pimpleKInventory=="WJNativeV1"); }
    bool ProspectiveKInventory() const { return PimpleKBoundEnabled() && m_pimpleKInventory=="ConstitutiveV3"; }
    std::string BranchCoverageDiagnostic(const Vector &state,double time,
        const RANSWallEvolutionData &data,const std::vector<RANSWallCoordinate> &anchors) const;
    std::string CheckConstitutiveCoverage(const Vector &state,double time,
        const RANSWallEvolutionData &data,const RANSWallCoordinate &unregistered) const;
    const std::string &PimpleKBoundMode() const { return m_pimpleKBound; }
    std::string PimpleKBoundSignature() const;
    std::string PimpleKBoundControls() const;
    static bool PimpleKBoundWorkParameter(const std::string &name);
    struct KBoundResult
    {
        Vector state;
        bool success=false,activated=false;
        int iterations=0,massFactorApplications=0,kktSolves=0;
        double integral=0;
        std::string report,failure;
    };
    std::string PimpleKRateBudget(const Vector &state,const Vector &history,double time,
        const RANSWallEvolutionData &data) const;
    std::string SourceBalanceDiagnostic(const Vector &state,double time,
        const RANSWallEvolutionData &data,const Vector &stations,const Vector &heights) const;
    std::string PimpleContactDiagnostic(const Vector &state,double time,
        const RANSWallEvolutionData &data,const Vector &steps,
        const std::vector<RANSWallCoordinate> &witnesses,bool exportMass,bool isolate) const;
    KBoundResult ApplyPimpleKBound(const Vector &raw,const Vector &preK,double time,
        const RANSWallEvolutionData &data,bool projectionOnly=false) const;
    bool SegregatedPicard() const { return m_nonlinearPolicy=="SegregatedPicard"; }
    Vector Norms(const Vector &residual,const std::map<int,double> &constraints,double ustar=.4) const;
    Vector StateNorms(const Vector &state,const Vector &reference,double ustar=.4) const;
    Vector ExplicitBody(const Vector &state,double time,const RANSWallEvolutionData &data) const;
    void SetReference(const Vector &reference,const RANSWallEvolutionData &data,double time,double eta);
    void RestoreReference(const Vector &reference,const Vector &defect,double eta,
        const RANSWallEvolutionData &data,double time);
    const Vector &Reference() const { return m_reference; }
    const Vector &Defect() const { return m_defect; }
    const Vector &RawReferenceLoads() const { return m_referenceLoads; }
    const std::map<int,double> &ReferenceConstraints() const { return m_referenceConstraints; }
    double Eta() const { return m_eta; }
    Vector DefectIntegrals() const;
    double AnchorEnergyWork(const Vector &state) const;
    std::vector<int> CondensableVelocity(const std::map<int,double> &constraints) const;
    double SampleK(const Vector &state,const RANSWallCoordinate &point) const;
    bool ColumnAverageK() const {return m_kSampling=="ColumnAverage";}
    bool HeightComparisonEnabled() const {return m_heightComparison || m_acceptedComparison;}
    bool AcceptedComparisonEnabled() const {return m_acceptedComparison;}
    Vector AcceptedPerturbation(const Vector &state,double amplitude) const;
    std::string FormulationMetadata() const;
    RANSWallGlobalPoint SampleKRow(const RANSWallCoordinate &anchor) const;
    // Exact physical mass restricted to the supplied global free-k rows.
    // Read-only access also supports native functional/projection verification.
    BlockMatrix PicardKMass(const std::vector<int> &rows) const;
    std::string HeightComparison(const Vector &state,const Vector &history,double time,
        const RANSWallEvolutionData &data,const Vector &reportingBreakpoints={}) const;
    std::string PhysicalKBands(const Vector &state,double time,const RANSWallEvolutionData &data,
        const Vector *correction=nullptr,bool exportIntegralRows=false) const;
    std::array<double,2> Forcing(const Vector &state,const RANSWallCoordinate &point,double time,
        const RANSWallEvolutionData &data,bool includeBody=true) const;
    std::array<double,3> CartesianForcing(const Vector &state,const RANSWallCoordinate &point,double time,
        const RANSWallEvolutionData &data,bool includeBody=true) const;
    Vector RoughnessBudget(const Vector &state,double time,const RANSWallEvolutionData &data) const;
    Vector EnergyBudget(const Vector &state,const Vector &history,double alpha,double time,const RANSWallEvolutionData &data) const;
    struct ResidualResolutionResult
    {
        Vector effectiveState;
        std::map<int,double> constraints;
        Diagnostic native,terms,compensated,absoluteAccumulation,additionCounts,norms;
        std::map<std::string,std::vector<long double>> extended;
        std::string scope;
    };
    ResidualResolutionResult ResidualResolution(const Vector &state,const Vector &history,
        double alpha,double time,const RANSWallEvolutionData &data,bool anchored=false) const;
    ResidualResolutionResult ResidualResolutionComponents(const Vector &state,const Vector &history,
        double alpha,double time,const RANSWallEvolutionData &data,bool anchored=false,
        bool originalOuterCounterfactual=false) const;
    /** Read-only constitutive inventory and one-sided tensor witnesses. No
     * lift, correction, governing solve or change to the selected model. */
    std::string ConstitutiveStateDiagnostic(const Vector &state,double time,
        const RANSWallEvolutionData &data,const std::vector<RANSWallCoordinate> &witnesses) const;
    std::string CompareOuterOperators(const Vector &state,double time,
        const RANSWallEvolutionData &data) const;
    /** Optional accepted-state audit. Entries are documented in the M7 guide.
     * Sampled admissibility is not a certificate over the entire polynomial. */
    Diagnostic Audit(const Vector &state,const Vector &history,double alpha,double time,
        const RANSWallEvolutionData &data,double ustar=.4,bool constitutiveReconstruction=false) const;
    const std::string &VariantName() const { return m_variantName; }
    bool CappedDestruction() const { return m_cappedDestruction; }
    bool SolveDiagnosticsEnabled() const { return m_solveDiagnostics; }
    bool SolveDiagnosticWritten() const { return m_solveDiagnosticWritten; }
    const std::string &SolveDiagnosticFile() const { return m_solveDiagnosticFile; }
    Vector RobustParameters() const;
    const SolverUtils::RANSVector &RobustReferenceFlow() const { return m_referenceFlow; }
    const SolverUtils::RANSVector &RobustReferenceNormal() const { return m_referenceNormal; }
    std::string InnerStressName() const;
    SolverUtils::RANSWallInnerStressValues InnerStress(double k,double auxiliaryNu,
        const SolverUtils::RANSTensor &gradient) const;
    SolverUtils::RANSClosureValues OuterClosure(double k,double epsilon,
        const SolverUtils::RANSTensor &gradient,const RANSWallCoordinate &point) const;
    SolverUtils::RANSCappedDestructionValues WallNormalization(const Vector &state,
        const RANSWallCoordinate &point,const SolverUtils::EquilibriumAWFProfiles &profile) const;
    SolverUtils::AWFProfileValues WallSource(const Vector &state,const RANSWallCoordinate &point,
        double production,const SolverUtils::EquilibriumAWFProfiles &profile,
        const SolverUtils::RANSCappedDestructionValues *normalization=nullptr,
        double *destructionRate=nullptr) const;
    double WallDissipation(const Vector &state,const RANSWallCoordinate &point,
        const SolverUtils::EquilibriumAWFProfiles &profile) const;
    std::string Signature() const;
    Diagnostic WJCompleteAudit(const Vector &state,const Vector &history,double alpha,double time,
        const RANSWallEvolutionData &data,bool rangesOnly=false) const;
    std::string WJAuditPacket(const Diagnostic &,const Vector &state,const Vector &history,
        double alpha,double time,const std::string &phase) const;
    std::string WJMeasurements(const Vector &state,const Vector &stageHistory,double stageAlpha,double time,
        const RANSWallEvolutionData &,const std::vector<std::array<double,8>> &observations,
        bool passiveDiagnostics=false) const;
    /** Passive pressure-test, divergence and signed near-roughness loads.
     * Any exported mass projection is offline diagnostic algebra, never a
     * governing correction or a modification of the accepted state. */
    std::string WJPassiveDiagnostics(const Vector &,const Vector &,double,double,
        const RANSWallEvolutionData &) const;
    std::string WJNativeChecks(const Vector &state,const RANSWallEvolutionData &data,double time,
        const Vector &acceptedHistory={},double acceptedAlpha=0) const;
    // Passive mapped checks. Frozen blocks are not a nonlinear Jacobian;
    // this operation performs no governing solve or initialization Lift.
    std::string WJTerrainAssemblyChecks(const Vector &state,const Vector &history,
        double alpha,double time,const RANSWallEvolutionData &data) const;
    std::string WJMappedBackendChecks(const Vector &,const Vector &,double,double,const RANSWallEvolutionData &) const;
    std::string WJWeakLoads(const Vector &,const Vector &,double,double,const RANSWallEvolutionData &) const;
    const std::string &ConstitutiveFailureContext() const { return m_solveInvalidContext; }
private:
    std::string m_boundaryMode="PrescribedV1",m_openTop="SlipV1",m_linearBackend="Dense";
    std::string m_mappedResourceDesign="None";
    std::string m_initialVelocityMode="Off";
    int m_initialVelocityCalls=0;
    class Factor
    {
    public:
        explicit Factor(const WallMatrix &matrix);
        explicit Factor(const BlockMatrix &matrix);
        Vector Apply(Vector rhs,bool transpose=false) const;
        std::string AllocationJSON() const;
    private:
        int n;
        Vector lu;
        std::vector<int> pivots;
        std::shared_ptr<RANSWallAlgebra::MappedBandFactor> band;
    };
    struct PicardBlock
    {
        std::vector<int> rows;
        BlockMatrix matrix;
        Vector residual;
    };
    void ExteriorPicard(int,const Vector &,double,const RANSWallEvolutionData &,
        PicardBlock &,const std::vector<int> &) const;
    struct PicardFlowResult
    {
        PicardBlock block;
        Vector predictor,step,defect;
        int velocitySize=0,pressureSize=0;
        double predictorDefect=0.;
    };
    PicardFlowResult PicardFlowCorrection(const Vector &state,const Vector &history,
        double alpha,double time,const RANSWallEvolutionData &data,
        const std::map<int,double> &external) const;
    Vector SolvePimpleSingleSweep(Vector state,const Vector &history,double alpha,double time,
        const RANSWallEvolutionData &data,int &iterations,bool denseReference=false) const;
    struct NonKResponse
    {
        bool converged=false;
        Vector state,residual,norms;
        int updates=0,flowSolves=0,epsilonSolves=0;
        double error=1e300;
        std::string failure;
    };
    NonKResponse EquilibrateNonK(const Vector &state,const Vector &history,double alpha,
        double time,const RANSWallEvolutionData &data,double tolerance,int maximumUpdates) const;
    Vector SolveReducedPicard(Vector state,const Vector &history,double alpha,double time,
        const RANSWallEvolutionData &data,int &iterations,double tolerance,int maximumIterations) const;
    Vector SolvePicardContinuation(Vector state,const Vector &history,double alpha,double time,
        const RANSWallEvolutionData &data,int &iterations,double tolerance,int maximumIterations) const;
    PicardBlock AssemblePicardBlock(int block,const Vector &state,const Vector &history,
        double alpha,double time,const RANSWallEvolutionData &data,
        const std::map<int,double> &external,bool denseReference=false) const;
    Vector ProjectPicardGround(const Vector &state,const PicardBlock &block,
        const Factor &factor,const Vector &step,double relaxation,int &cuts,
        std::string *diagnostic=nullptr) const;
    Vector SolveSegregatedPicard(Vector state,const Vector &history,double alpha,double time,
        const RANSWallEvolutionData &data,int &iterations,double tolerance,int maximumIterations) const;
    struct WallExtremum
    {
        RANSWallGlobalPoint row;
        long double value,error;
    };
    std::vector<WallExtremum> WallTraceR4(const Vector &state) const;
    std::string WallTraceJSON(const Vector &state,bool &admissible) const;
    bool RestrictedInner(const ScalarObservation &p) const;
    RANSWallGlobalPoint RestrictedGradientRow(const ScalarObservation &p) const;
    struct OuterBound
    {
        RANSWallGlobalPoint row;
        int family,block,axis,side,observationGeometry;
        bool projected;
        double scale,value,reference,width,lower,upper;
        std::string inventory;
        double Action(const Vector &step,int offset) const;
    };
    struct OuterViolation
    {
        std::array<double,4> family{};
        double sum=0;
        double Max() const;
        std::string report;
    };
    std::vector<OuterBound> OuterBoundsR5(const ScalarDomainResult &domain) const;
    OuterViolation OuterViolationR5(const std::vector<OuterBound> &bounds,
        const Vector &effective={},const Vector &trace={},double fraction=0.) const;
    Vector CoupledDirectionR4(const ScalarDomainResult &domain,const Vector &newton,
        const WallMatrix &jac,const std::map<int,double> &constraints,
        const std::function<Vector(Vector,bool)> &inverse,std::string &details,const std::vector<OuterBound> &outer={},
        double restoration=1.,bool inward=true,
        RANSWallR3::DurableProgress *diagnostic=nullptr,const std::string &tag="",const Vector &residual={},
        double radius=2.,double *guardWork=nullptr,const std::array<double,4> &interior={}) const;
    WallMatrix AssembleNonpressureR3(const Vector &state,double time,
        const RANSWallEvolutionData &data,bool imposeConstraints,int region,
        std::vector<std::string> *metadata,const std::vector<int> *freeK=nullptr) const;
    WallMatrix PicardKJacobian(const Vector &state,double alpha,double time,
        const RANSWallEvolutionData &data,const std::vector<int> &rows) const;
    Vector PicardKResidual(const Vector &state,const Vector &history,double alpha,double time,
        const RANSWallEvolutionData &data,const std::vector<int> &rows) const;
    void DiagnosePicardK(const Vector &state,const Vector &history,double alpha,double time,
        const RANSWallEvolutionData &data,const std::map<int,double> &external) const;
    struct ScalarLocation
    {
        std::string inventory;
        RANSWallGlobalPoint row;
        int block=0,side=0;
        bool strict=false;
    };
    void InitialiseScalarLocations() const;
    RANSWallSharedInventory<ScalarLocation> ReconstructionLocations(int block,bool complete) const;
    std::vector<RANSWallCoordinate> ConstitutiveAnchors(const Vector &state,
        const std::vector<RANSWallCoordinate> &additional) const;
    std::vector<std::pair<RANSWallCoordinate,double>> PhysicalBandRays() const;
    std::vector<RANSWallGlobalPoint> NormalizationRows(const RANSWallCoordinate &anchor,
        const SolverUtils::EquilibriumAWFProfiles &profile) const;
    std::shared_ptr<SolverUtils::EquilibriumAWFProfiles> ProspectiveActiveProfile(
        const Vector &state,const RANSWallCoordinate &anchor,std::string &reason) const;
    // A coverage tripwire in the final audit checks actual consumers, in
    // addition to sharing their location generators with the projection.
    struct CoverageScope
    {
        bool &flag,previous;
        CoverageScope(bool &value,bool enable):flag(value),previous(value){flag=enable;}
        ~CoverageScope(){flag=previous;}
    };
    mutable bool m_assertKCoverage=false;
    mutable std::set<RANSWallCoordinate> m_coveredKLocations;
    mutable std::set<RANSWallCoordinate> m_coveredKFunctionals;
    std::string m_kSampling="Point";
    bool m_heightComparison=false,m_acceptedComparison=false;
    mutable std::map<RANSWallCoordinate,RANSWallGlobalPoint> m_kSamplingRows;
    bool m_evaluationReuse=false;
    mutable std::shared_ptr<EvaluationCache> m_evaluation;
    mutable std::map<std::pair<int,bool>,RANSWallSharedInventory<ScalarLocation>> m_reconstructionCatalog;
    mutable std::map<std::tuple<int,int,double,int,int>,std::vector<FacePoint>> m_surfaceCatalog;
    mutable std::vector<std::map<int,int>> m_geometryCatalog;
    struct EvaluationWork
    {
        std::uint64_t scopes=0,domainCalls=0,domainComputes=0,domainHits=0,prefixHits=0,
            rowConstructs=0,rowHits=0,rowEvictions=0,reconstructionBuilds=0,reconstructionHits=0,
            surfaceBuilds=0,surfaceHits=0,normalizationCalls=0,normalizationComputes=0,normalizationHits=0,
            residualCalls=0,residualHits=0,constraintCalls=0,constraintHits=0,
            stateEvictions=0,inputInvalidations=0,profileReuses=0,scalarFieldReuses=0,peakStates=0,peakRows=0,peakProfiles=0;
    };
    mutable EvaluationWork m_evaluationWork;
    void PrepareEvaluation(const RANSWallEvolutionData &,double) const;
    std::shared_ptr<const RANSWallGlobalPoint> EvaluationRow(int block,int element,const RANSWallCoordinate &) const;
    int GeometryElement(int block,int geometry) const;
    std::shared_ptr<const ScalarDomainResult> CachedDomain(const Vector &,double,const RANSWallEvolutionData &,bool,
        const std::vector<RANSWallCoordinate> &) const;
    void StoreDomain(const Vector &,double,const RANSWallEvolutionData &,bool,const std::vector<RANSWallCoordinate> &,
        const ScalarDomainResult &) const;
    void InvalidateCoverageToken() const;
    ScalarDomainResult EvaluateScalarDomain(const Vector &,double,const RANSWallEvolutionData &,bool,bool,
        const std::vector<RANSWallCoordinate> &) const;
    std::shared_ptr<SolverUtils::EquilibriumAWFProfiles> EvaluationProfile(
        const std::shared_ptr<SolverUtils::EquilibriumAWFProfiles> &,double) const;
    std::shared_ptr<const SolverUtils::RANSCappedDestructionValues> CachedNormalization(const Vector &,
        const RANSWallCoordinate &,const SolverUtils::EquilibriumAWFProfiles &) const;
    void StoreNormalization(const Vector &,const RANSWallCoordinate &,const SolverUtils::EquilibriumAWFProfiles &,
        const SolverUtils::RANSCappedDestructionValues &) const;
    std::shared_ptr<const Vector> CachedResidual(const Vector &,const Vector &,double,double,
        const RANSWallEvolutionData &,bool,bool,int) const;
    void StoreResidual(const Vector &,const Vector &,double,double,const RANSWallEvolutionData &,bool,bool,int,const Vector &) const;
    std::shared_ptr<const std::map<int,double>> CachedConstraints(const Vector &,double,const RANSWallEvolutionData &) const;
    void StoreConstraints(const Vector &,double,const RANSWallEvolutionData &,const std::map<int,double> &) const;
    Vector SolveScalarFeasible(Vector state,const Vector &history,double alpha,
        double time,const RANSWallEvolutionData &data,bool condensed,int &iterations,
        double tolerance,int maximumIterations,bool scalarOnly=false) const;
    struct VolumePoint { RANSWallGlobalPoint v,p,e; double weight; };
    struct Ray
    {
        RANSWallCoordinate coordinate;double weight,cachedK=-1,cachedWeight=-1;
        RANSWallGlobalPoint sample,outerV,outerP,outerK;
        std::shared_ptr<SolverUtils::EquilibriumAWFProfiles> profile;
        SolverUtils::RANSCappedDestructionValues destruction;
        std::vector<VolumePoint> points;
    };
    std::map<int,double> ExternalConstraints(double time,const RANSWallEvolutionData &data) const;
    Vector OwnedBodyLoad(double time,const RANSWallEvolutionData &data,int region=-1) const;
    std::shared_ptr<SolverUtils::EquilibriumAWFProfiles> Profile(double k,const RANSWallCoordinate &chart={},double weight=-1) const;
    RANSWallGlobalPoint Transform(RANSWallGlobalPoint point) const;
    BlockMatrix PhysicalMass(int block,int region=-1) const;
    Vector PhysicalProject(int block,const RANSWallEvolutionFunction &function,double time) const;
    void Refresh(const Vector &state,double time=0,const RANSWallEvolutionData *data=nullptr) const;
    Vector ExternalLift(const Vector &state,double time,const RANSWallEvolutionData &data) const;
    void InitialiseTransition(const LibUtilities::SessionReaderSharedPtr &session);
    void InitialiseRobust(const LibUtilities::SessionReaderSharedPtr &session);
    Vector FlatReference(double ustar,const std::array<double,2> &direction) const;
    void Broadcast(Vector &value) const;
    struct SolverCache;
    int m_dimension;
    LibUtilities::CommSharedPtr m_comm;
    RANSWallPatchParameters m_parameters;
    SolverUtils::RANSWallNormalMap m_mapping;
    int m_ownershipShift=0;
    SolverUtils::AWFProfileOptions m_options;
    std::string m_formulation="EPAWF";
    std::string m_nativeGeometry="FlatV1";
    SolverUtils::AWFParameters m_wallParameters;
    std::string m_variantName="Original";
    SolverUtils::RANSRobustVariant m_robust;
    SolverUtils::RANSCappedDestructionControls m_destruction;
    bool m_cappedDestruction=false;
    SolverUtils::RANSWallInnerStressModel m_innerStress=SolverUtils::RANSWallInnerStressModel::Linear;
    double m_regularizationGamma=0,m_regularizationScale=0,m_robustUstar=.4;
    SolverUtils::RANSVector m_referenceFlow{1,0,0},m_referenceNormal{0,0,1};
    bool m_transitionEnabled=false,m_transitionReady=false,m_uniformEP=false;
    SolverUtils::AWFTransitionControls m_transition;
    double m_sensorA=2.5,m_sensorB=3.5,m_baselineStart=-1.25;
    Vector m_transitionReference;
    std::set<int> m_traceIds,m_baselineTraceIds;
    using ProfilePair=std::array<std::shared_ptr<SolverUtils::EquilibriumAWFProfiles>,2>;
    mutable std::map<std::array<double,3>,ProfilePair> m_profileCache;
    mutable std::map<std::array<double,3>,std::shared_ptr<SolverUtils::EquilibriumAWFProfiles>> m_prospectiveProfileCache;
    std::vector<BlockMatrix> m_mass;
    mutable std::array<BlockMatrix,2> m_regionalMass;
    mutable Vector m_boundaryFlux,m_boundaryDirection;
    mutable double m_boundaryCost=0;
    mutable double m_boundaryCorrection=0;
    Vector m_volumes;
    std::vector<Vector> m_constants;
    std::array<Vector,2> m_pressureMean;
    mutable std::vector<Ray> m_rays;
    std::vector<VolumePoint> m_outer,m_bodyInner;
    std::vector<int> m_boundaryTopology;
    Vector m_reference,m_defect,m_referenceLoads;
    std::map<int,double> m_referenceConstraints;
    double m_eta=0;
    mutable std::shared_ptr<SolverCache> m_solver;
    // Opt-in, serial first-failure diagnostics. A counterfactual same-stage
    // fresh-Jacobian retry is recorded separately and never advances history.
    bool m_solveDiagnostics=false,m_solveDiagnosticReplay=false,m_solveDiagnosticAdmissible=false;
    int m_solveDiagnosticIterations=12;
    std::string m_solveDiagnosticFile;
    mutable bool m_solveDiagnosticWritten=false;
    mutable std::string m_solveResidualException;
    mutable std::string m_solveInvalidContext;
    std::string m_nonlinearPolicy="Legacy",m_feasibleTraceFile="case.feasible.jsonl";
    std::string m_temporalCoupling="ConvergedStage";
    std::string m_denseCapacityRequest="3000",m_denseResourceDesign="None";
    int m_denseCapacity=3000;
    bool m_pimpleVerify=false;
    std::string m_pimpleKBound="Off",m_pimpleSnapshotDirectory="pimple-states";
    std::string m_pimpleKInventory="LegacyV1";
    double m_pimpleKFloorRelative=1e-8,m_pimpleBudgetTime=-1.;
    bool m_pimpleSnapshots=false,m_pimpleArithmetic=false;
    double m_pimpleArithmeticTime=-1;
    int m_kBoundMaxAdd=2000,m_kBoundMaxTotal=2000,m_kBoundMaxSeparations=2000;
    double m_kBoundMaxSeconds=900.;
    std::string m_kBoundTrace="Summary",m_kBoundTraceFile="None",m_kBoundDumpFile="None";
    std::string m_kBoundSolver="WorkingSet";
    mutable int m_liftCalls=0,m_liftIterations=0;
    mutable std::vector<std::string> m_pimpleRecords;
    std::array<double,3> m_picardRelaxation{{1.,.5,.5}};
    double m_picardPseudoAlpha=0.;
    double m_picardGuessPerturbation=0.;
    double m_picardNonKTolerance=0.;
    bool m_picardVerify=false;
    bool m_picardGroundProjection=false;
    bool m_picardKNewton=false;
    bool m_picardReduced=false,m_reducedVerify=false;
    bool m_picardContinuation=false;
    double m_picardContinuationEta=0.,m_picardContinuationTolerance=0.;
    int m_reducedBasis=10,m_reducedHUpdates=50;
    double m_reducedRadius=.05,m_reducedStep=2e-5,m_reducedHTolerance=2e-9;
    bool m_picardEquilibrateK=false;
    double m_picardKRadius=.15;
    std::string m_picardGuessFile,m_picardKDiagnostic;
    std::string m_directionAccuracy="Legacy",m_r7Policy="Off";
    double m_r7PreviousDisplacement=0.,m_r7Radius=0.;
    bool m_r7Interior=false;
    std::string m_r3ProgressFile="case.r3-progress.jsonl";
    std::string m_actionPolicy="Legacy";
    std::string m_globalization="FractionToBoundary";
    mutable bool m_feasibleTraceStarted=false;
    mutable std::vector<ScalarLocation> m_scalarLocations;
    mutable std::vector<RANSWallCoordinate> m_traceScalarAnchors,m_scalarSamplerAnchors;
    mutable Vector m_lastFeasibleState;
    mutable int m_feasibleSolveCount=0;
};
}
#endif
