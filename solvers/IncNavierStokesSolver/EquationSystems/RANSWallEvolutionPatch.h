///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_RANS_WALL_EVOLUTION_PATCH_H
#define NEKTAR_RANS_WALL_EVOLUTION_PATCH_H
#include <IncNavierStokesSolver/EquationSystems/RANSWallGlobalField.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallPatch.h>
#include <LibUtilities/Communication/Comm.h>
#include <SolverUtils/Turbulence/RANSWallNormalMap.h>
#include <SolverUtils/Turbulence/EquilibriumAWFTransition.h>
#include <SolverUtils/Turbulence/RANSRobustVariant.h>
#include <SolverUtils/Turbulence/RANSWallInnerStress.h>
#include <set>

namespace Nektar
{
namespace RANSWallR3 { class DurableProgress; }
SolverUtils::RANSWallNormalMap LoadRANSWallNormalMap(const LibUtilities::SessionReaderSharedPtr &session,int dimension,double height);
using RANSWallEvolutionFunction=std::function<double(const RANSWallCoordinate&,double)>;
struct RANSWallEvolutionData
{
    std::vector<RANSWallEvolutionFunction> innerBody,outerBody,boundary;
    // Sensors retain the physical body force when an IMEX stage moves its
    // prescribed load to the explicit side. Empty means use outerBody.
    std::vector<RANSWallEvolutionFunction> sensorBody;
    RANSWallEvolutionFunction innerKSource,outerKSource,outerEpsilonSource;
    // Optional fixed native external traces, e.g. the stored neutral target.
    std::map<int,double> fixedBoundary;
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
public:
    using Vector=RANSWallAlgebra::Vector;
    using WallMatrix=RANSWallAlgebra::DenseMatrix;
    using Diagnostic=std::map<std::string,Vector>;
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
    int Size() const { return offset[Gauge()]+1; }
    int PhysicalAxis(int a) const { return m_dimension==2 && a==1?2:a; }
    const RANSWallPatchParameters &Parameters() const { return m_parameters; }
    const LibUtilities::CommSharedPtr &Comm() const { return m_comm; }
    const SolverUtils::RANSWallNormalMap &Mapping() const { return m_mapping; }
    bool Mapped() const { return !m_mapping.Identity(); }
    int Owner(int geometry) const { return (geometry+m_ownershipShift)%m_comm->GetSize(); }
    RANSWallCoordinate Physical(const RANSWallCoordinate &chart) const;
    RANSWallGlobalPoint Point(int block,const RANSWallCoordinate &chart,int side=0) const;
    std::vector<RANSWallGlobalPoint> Quadrature(int block,int order=0) const;
    struct FacePoint { RANSWallGlobalPoint test;SolverUtils::RANSVector area;double weight; };
    std::vector<FacePoint> Surface(int block,int axis,double position,int side,int order=0) const;
    double Distance(const RANSWallCoordinate &chart) const;
    double DependentEpsilon(const Vector &state,const RANSWallCoordinate &chart) const;
    double DependentEpsilon(const Vector &state,const RANSWallCoordinate &chart,double time,
        const RANSWallEvolutionData &data) const;
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
    };
    struct ScalarDomainResult
    {
        bool admissible=false;
        std::string phase,failure,limitingInventory,profileBranches;
        RANSWallCoordinate limitingPoint{};
        double minimumMargin=1e300,minimumInnerK=1e300,minimumOuterK=1e300,
            minimumOuterEpsilon=1e300,minimumSampledK=1e300,minimumTraceEpsilon=1e300,
            minimumB=1,minimumLambda=1,minimumCutoff=1e300,maximumCutoff=-1e300;
        std::map<std::string,int> counts;
        Vector effectiveState,traceState;
        std::vector<ScalarObservation> observations;
    };
    // Scalar evaluation domain only: covariance and the restricted outer
    // physical box are separate accepted-stage gates.
    ScalarDomainResult ScalarDomain(const Vector &state,double time,
        const RANSWallEvolutionData &data,bool retainLocations=false) const;
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
    const std::string &TemporalCouplingName() const { return m_temporalCoupling; }
    Vector InitializationWork() const { return {double(m_liftCalls),double(m_liftIterations)}; }
    const std::vector<std::string> &PimpleRecords() const { return m_pimpleRecords; }
    Diagnostic PhysicalSamples(const Vector &state,double time,const RANSWallEvolutionData &data) const;
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
        double alpha,double time,const RANSWallEvolutionData &data,bool anchored=false) const;
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
private:
    class Factor
    {
    public:
        explicit Factor(const WallMatrix &matrix);
        Vector Apply(Vector rhs,bool transpose=false) const;
    private:
        int n;
        Vector lu;
        std::vector<int> pivots;
    };
    struct PicardBlock
    {
        std::vector<int> rows;
        WallMatrix matrix;
        Vector residual;
    };
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
        const RANSWallEvolutionData &data,int &iterations) const;
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
        const std::map<int,double> &external) const;
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
    WallMatrix PicardKMass(const std::vector<int> &rows) const;
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
    WallMatrix PhysicalMass(int block,int region=-1) const;
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
    std::string m_variantName="Original";
    SolverUtils::RANSRobustVariant m_robust;
    SolverUtils::RANSCappedDestructionControls m_destruction;
    bool m_cappedDestruction=false;
    SolverUtils::RANSWallInnerStressModel m_innerStress=SolverUtils::RANSWallInnerStressModel::Linear;
    double m_regularizationGamma=0,m_regularizationScale=0,m_robustUstar=.4;
    SolverUtils::RANSVector m_referenceFlow{1,0,0},m_referenceNormal{0,0,1};
    bool m_transitionEnabled=false,m_transitionReady=false;
    SolverUtils::AWFTransitionControls m_transition;
    double m_sensorA=2.5,m_sensorB=3.5,m_baselineStart=-1.25;
    Vector m_transitionReference;
    std::set<int> m_traceIds,m_baselineTraceIds;
    using ProfilePair=std::array<std::shared_ptr<SolverUtils::EquilibriumAWFProfiles>,2>;
    mutable std::map<std::array<double,3>,ProfilePair> m_profileCache;
    std::vector<WallMatrix> m_mass;
    mutable std::array<WallMatrix,2> m_regionalMass;
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
    bool m_pimpleVerify=false;
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
