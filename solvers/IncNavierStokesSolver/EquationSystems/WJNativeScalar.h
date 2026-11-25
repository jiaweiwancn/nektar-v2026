// SPDX-License-Identifier: MIT
#ifndef NEKTAR_WJ_NATIVE_SCALAR_H
#define NEKTAR_WJ_NATIVE_SCALAR_H

#include <IncNavierStokesSolver/EquationSystems/WJNativeRules.h>
#include <array>
#include <unordered_map>

namespace Nektar::WJNative
{
using ScalarBoundaryFunction =
    std::function<double(double, double, const Vector &)>;

struct ScalarSolveInfo
{
    int iterations = 0, rebuilds = 0;
    double residual = 0;
    double factorSeconds = 0, solveSeconds = 0, residualSeconds = 0;
};

/** Native continuous Q3 scalar algebra on the full flow domain or a rectangular
 * geometry subset. The scalar signed assembly map is independent of u, v and p.
 * This class supplies no turbulence model, clipping, scalar stabilisation or
 * physical coefficient update. All diffusion/reaction/source samples are inputs.
 */
class Scalar
{
public:
    /** A nonempty essentialSides list defines the complete essential mask:
     * 0 left, 1 right, 2 bottom, 3 top, relative to THIS field's rectangle.
     * In particular epsilon's bottom is the interface, not the ground.
     * An empty list uses the field's native Dirichlet coefficient mask.
     * residualScale must be positive and is a physical reference independent
     * of the time step and the current scalar/linear-system diagonal.
     */
    Scalar(const MultiRegions::ContFieldSharedPtr &field, const Flow &flow,
           Function residualScale, std::vector<int> essentialSides = {});

    int Size() const { return m_size; }
    int DirichletSize() const { return int(m_essentialIds.size()); }
    int FreeSize() const { return int(m_freeIds.size()); }
    const std::vector<int> &FreeIds() const { return m_freeIds; }
    const MultiRegions::ContFieldSharedPtr &Field() const { return m_field; }
    const Sparse &Mass() const { return m_mass; }
    bool Owns(int geometry) const { return m_geometryToCell.count(geometry)!=0; }
    const std::vector<int> &CellIds(int geometry) const
    { return m_cells.at(m_geometryToCell.at(geometry)).ids; }
    Vec Evaluate(const NativeRule &rule,const Vec &state,int derivative=-1) const;
    /** Sufficient full-cell bounds from the native Q3 tensor polynomial's
     * Bernstein coefficients, including endpoints and inter-node regions. */
    std::array<double,2> PolynomialBounds(const Vec&state,int subdivisionDepth=0)const;
    void ResetPreconditioner() { m_factor.reset(); m_initialIterations=m_lastIterations=0; }

    Vec Project(const Function &function) const;
    /** Endpoints are interpolated exactly; only the two interior edge modes
     * are L2 projected. Repeated essential corner values must agree.
     */
    Vec EssentialLift(const ScalarBoundaryFunction &function,
                      const RuleFactory &factory) const;
    Vec Lift(const Vec &fullDirValues) const;
    Vec Restrict(const Vec &full) const;
    Vec Expand(const Vec &free, const Vec &lift) const;

    Sparse System(double alpha, const Sparse &K) const;
    Vec LiftedRhs(double alpha, const Sparse &K, const Vec &fullLoad,
                  const Vec &lift) const;
    Sparse Diffusion(const Rules &rules, const std::vector<Vec> &diffusion,
                     const std::vector<Vec> &reaction = {}) const;
    /** Full-row current diffusion/reaction action with this scalar's signed
     * assembly map. It equals Diffusion(rules,D,reaction)*state without
     * constructing the sparse matrix; essential reactions are retained. */
    Vec DiffusionAction(const Rules &rules, const std::vector<Vec> &diffusion,
                        const Vec &state,
                        const std::vector<Vec> &reaction = {}) const;
    Vec Source(const Rules &rules, const std::vector<Vec> &source) const;
    /** Full-row conservative LHS load: -int grad(phi).U q + int_boundary phi Un q.
     * Velocity uses canonical [u,v] coefficients (an optional p tail is ignored).
     */
    Vec Advection(const Rules &rules, const Vec &fullCanonicalVelocity,
                  const Vec &ownGlobalScalar) const;
    /** Full native conservative-transport matrix with a supplied frozen
     * velocity. This is also useful as a steady preconditioner block. */
    Sparse TransportMatrix(const Rules &rules,const Vec &fullCanonicalVelocity)const;
    /** RHS load -int_boundary phi j_out, where j_out=-D grad(q).n.
     * The callback is visited on every exterior face, including essential
     * faces. The caller chooses physical fluxes; this routine invents none.
     */
    Vec NaturalFlux(const ScalarBoundaryFunction &outwardFlux,
                    const RuleFactory &factory,
                    const std::function<double(const Element&,double)>&cutoff={}) const;

    ScalarSolveInfo Solve(const Sparse &current, const Vec &rhs, Vec &free,
                          double coefficientChange = 0,
                          double tolerance = 1e-11);
    ScalarSolveInfo Residual(const Sparse &current, const Vec &rhs,
                             const Vec &free) const;
    /** RMS dual norm over the homogeneous test space:
     * sqrt(r_f^T [int scale^2 phi_f phi_f^T]^{-1} r_f / volume).
     * Essential-row reactions are omitted. There is no alpha normalisation.
     */
    double Dual(const Vec &fullResidual) const;
    double FreeResidualNorm(const Vec&residual)const{return FreeDual(residual);}
    const Vec&ResidualRowScale()const{return m_scale;}
    double PrepareLinearPreconditioner(const Sparse&current){return Prepare(current);}
    Vec ApplyLinearPreconditioner(const Vec&rhs)const{return Precondition(rhs);}

private:
    struct Cell
    {
        const Element *flow = nullptr;
        std::vector<int> ids;
        std::vector<double> signs;
    };
    struct Face
    {
        int cell = 0, side = 0;
        NativeRule rule;
    };

    const Cell *Find(const NativeRule &rule) const;
    Mat Basis(const NativeRule &rule, const Cell &cell, int derivative,
              bool integration=true) const;
    Vec Local(const Vec &full, const std::vector<int> &ids) const;
    void Add(Vec &full, const Vec &local, const Cell &cell) const;
    Sparse Gather(const Sparse &full) const;
    double FreeDual(const Vec &residual) const;
    double Prepare(const Sparse &current);
    Vec Precondition(const Vec &rhs) const;

    MultiRegions::ContFieldSharedPtr m_field;
    const Flow &m_flow;
    Function m_referenceScale;
    int m_size = 0;
    double m_volume = 0;
    bool m_explicitEssential = false;
    std::array<bool, 4> m_essentialSides{};
    std::vector<Cell> m_cells;
    std::vector<Face> m_faces;
    std::unordered_map<int, int> m_geometryToCell;
    std::vector<int> m_essentialIds, m_freeIds, m_fullToFree;
    Sparse m_mass;
    Vec m_scale;
    Eigen::SimplicialLDLT<Sparse> m_riesz;
    std::unique_ptr<Eigen::SparseLU<Sparse>> m_factor;
    int m_initialIterations = 0, m_lastIterations = 0;
};
}
#endif
