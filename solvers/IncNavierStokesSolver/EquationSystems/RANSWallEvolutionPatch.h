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
#include <set>

namespace Nektar
{
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
    Vector DirectionalAction(const Vector &state,const Vector &direction,const Vector &history,
        double alpha,double time,const RANSWallEvolutionData &data,double step=2e-6) const;
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
    /** Optional accepted-state audit. Entries are documented in the M7 guide.
     * Sampled admissibility is not a certificate over the entire polynomial. */
    using Diagnostic=std::map<std::string,Vector>;
    Diagnostic Audit(const Vector &state,const Vector &history,double alpha,double time,
        const RANSWallEvolutionData &data,double ustar=.4) const;
    std::string Signature() const;
private:
    struct VolumePoint { RANSWallGlobalPoint v,p,e; double weight; };
    struct Ray
    {
        RANSWallCoordinate coordinate;double weight,cachedK=-1,cachedWeight=-1;
        RANSWallGlobalPoint sample,outerV,outerP,outerK;
        std::shared_ptr<SolverUtils::EquilibriumAWFProfiles> profile;
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
    Vector FlatReference(double ustar,const std::array<double,2> &direction) const;
    void Broadcast(Vector &value) const;
    struct SolverCache;
    int m_dimension;
    LibUtilities::CommSharedPtr m_comm;
    RANSWallPatchParameters m_parameters;
    SolverUtils::RANSWallNormalMap m_mapping;
    int m_ownershipShift=0;
    SolverUtils::AWFProfileOptions m_options;
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
};
}
#endif
