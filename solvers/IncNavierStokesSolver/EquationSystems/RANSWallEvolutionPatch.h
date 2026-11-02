///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_RANS_WALL_EVOLUTION_PATCH_H
#define NEKTAR_RANS_WALL_EVOLUTION_PATCH_H
#include <IncNavierStokesSolver/EquationSystems/RANSWallGlobalField.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallPatch.h>
#include <LibUtilities/Communication/Comm.h>

namespace Nektar
{
using RANSWallEvolutionFunction=std::function<double(const RANSWallCoordinate&,double)>;
struct RANSWallEvolutionData
{
    std::vector<RANSWallEvolutionFunction> innerBody,outerBody,boundary;
    RANSWallEvolutionFunction innerKSource,outerKSource,outerEpsilonSource;
    // Optional fixed native external traces, e.g. the stored neutral target.
    std::map<int,double> fixedBoundary;
    explicit RANSWallEvolutionData(int dimension);
};

/** Flat 2D/3D native CG stage with globally ordered replicated coefficients.
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
        const LibUtilities::CommSharedPtr &collective,const RANSWallPatchParameters &parameters);
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
        const RANSWallEvolutionData &data,bool anchored=true,bool imposeConstraints=true) const;
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
        const RANSWallEvolutionData &data) const;
    std::string Signature() const;
private:
    struct VolumePoint { RANSWallGlobalPoint v,p,e; double weight; };
    struct Ray
    {
        RANSWallCoordinate coordinate;double weight,cachedK=-1;
        RANSWallGlobalPoint sample,outerV,outerP,outerK;
        std::shared_ptr<SolverUtils::EquilibriumAWFProfiles> profile;
        std::vector<VolumePoint> points;
    };
    std::map<int,double> ExternalConstraints(double time,const RANSWallEvolutionData &data) const;
    std::shared_ptr<SolverUtils::EquilibriumAWFProfiles> Profile(double k) const;
    void Refresh(const Vector &state) const;
    void Broadcast(Vector &value) const;
    struct SolverCache;
    int m_dimension;
    LibUtilities::CommSharedPtr m_comm;
    RANSWallPatchParameters m_parameters;
    SolverUtils::AWFProfileOptions m_options;
    std::vector<WallMatrix> m_mass;
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
