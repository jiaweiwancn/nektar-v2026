///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_RANS_WALL_PATCH_H
#define NEKTAR_RANS_WALL_PATCH_H
#include <IncNavierStokesSolver/EquationSystems/RANSWallFieldLayout.h>
#include <SolverUtils/Turbulence/EquilibriumAWFProfiles.h>
#include <memory>

namespace Nektar
{
using RANSWallTimeFunction = std::function<double(double,double,double)>;
struct RANSWallPatchData
{
    // Prescribed physical body force, sampled from the OUTER side in F_A.
    std::array<RANSWallTimeFunction,2> innerBody, outerBody;
    RANSWallTimeFunction innerKSource, outerKSource, outerEpsilonSource;
    // U,W,k,epsilon on external Dirichlet boundaries. Ground velocity is zero;
    // ground k flux is zero. Epsilon at H is supplied by the dependent profile.
    std::array<RANSWallTimeFunction,4> boundary;
    RANSWallPatchData();
};
struct RANSWallPatchParameters
{
    SolverUtils::AWFProfileGeometry geometry{1.5e-5,2,.1,1,1};
    std::string modelCalibration,modelManifest; // empty preserves historical CLS identity
    double chi = 1;
    int tangentialOrder = 8, normalOrder = 16;
    SolverUtils::RANSNeutralReference model = SolverUtils::RANSNeutralReference::CLS();
};

/** Uncondensed native 2D flat patch: conservative full-strain momentum and
 * scalar transport, mixed pressure with one global gauge, original CLS outside
 * and the revised EP-AWF inside. alpha*(state-history) retains inner storage.
 * This stage operator is independent of the later native IMEX time-map gate.
 */
class RANSWallPatch
{
public:
    using Vector = RANSWallAlgebra::Vector;
    using WallMatrix = RANSWallAlgebra::DenseMatrix;
    using Layout = RANSWallFieldLayout;
    RANSWallPatch(std::shared_ptr<Layout> layout, const RANSWallPatchParameters &parameters,
                 bool prepareMixedStage=true);
    const std::shared_ptr<Layout> &Fields() const { return m_layout; }
    const RANSWallPatchParameters &Parameters() const { return m_parameters; }
    Vector Project(const std::array<RANSWallTimeFunction,6> &functions, double time) const;
    double Forcing(const Vector &state, double x, double time, const RANSWallPatchData &data) const;
    double SampleK(const Vector &state, double x) const;
    std::map<int,double> Constraints(const Vector &state, double time, const RANSWallPatchData &data) const;
    Vector Residual(const Vector &state, const Vector &history, double alpha,
                    double time, const RANSWallPatchData &data) const;
    WallMatrix Jacobian(const Vector &state, const Vector &history, double alpha,
                        double time, const RANSWallPatchData &data, double step = 2e-6) const;
    Vector Solve(Vector state, const Vector &history, double alpha, double time,
                 const RANSWallPatchData &data, bool condensed, int &iterations,
                 double tolerance = 2e-12, int maximumIterations = 15) const;
    Vector ResidualNorms(const Vector &residual, const std::map<int,double> &constraints) const;
    std::vector<int> CondensableVelocity(const std::map<int,double> &constraints) const;
    const std::array<Vector,2> &PressureMeanLoads() const { return m_pressureMean; }
    void Store(const Vector &state) const;
    struct ParallelResult
    {
        Vector coefficients,metrics;
        WallMatrix samples;
        int radialUnknowns=0;
    };
    // Exact parallel/quasi-steady reduction only, with K and U_H fixed.
    // Native 2D CG functions are restricted to constant tangential modes.
    ParallelResult ParallelReduction(double sampledK,double matchingVelocity,double forcing) const;
private:
    struct Point { RANSWallPoint v,k,p,e; double weight; };
    struct Ray { double x,weight,k = -1; std::shared_ptr<SolverUtils::EquilibriumAWFProfiles> profile;
                 std::vector<Point> points; };
    std::shared_ptr<SolverUtils::EquilibriumAWFProfiles> Profile(double k) const;
    std::map<int,double> FixedConstraints(double time,const RANSWallPatchData &data) const;
    Vector Lift(const Vector &state,double time,const RANSWallPatchData &data,
                std::map<int,double> &constraints) const;
    std::shared_ptr<Layout> m_layout;
    RANSWallPatchParameters m_parameters;
    SolverUtils::AWFProfileOptions m_options;
    mutable std::vector<Ray> m_rays;
    std::vector<Point> m_outer;
    std::array<Vector,2> m_pressureMean;
    std::array<WallMatrix,6> m_mass;
    bool m_stagePrepared;
};

/** Exact algebraic Schur reduction of a verified uncondensed Newton matrix.
 * Only free inner velocity modes are removed. All pressure modes, including
 * the relative subdomain pressure mode, remain until traction matching/gauge.
 */
class RANSWallCondensation
{
public:
    using Vector = RANSWallAlgebra::Vector;
    using WallMatrix = RANSWallAlgebra::DenseMatrix;
    RANSWallCondensation(const WallMatrix &jacobian,const Vector &residual,const std::vector<int> &inner);
    Vector Recover(const Vector &retainedStep) const;
    Vector Solve() const;
    const WallMatrix &Matrix() const { return m_schur; }
    const Vector &Residual() const { return m_residual; }
    const std::vector<int> &Retained() const { return m_retained; }
private:
    std::vector<int> m_inner,m_retained;
    WallMatrix m_ii,m_ir,m_schur;
    Vector m_innerResidual,m_residual;
};
}
#endif
