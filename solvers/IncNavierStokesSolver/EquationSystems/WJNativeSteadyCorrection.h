// SPDX-License-Identifier: MIT
#ifndef NEKTAR_WJ_NATIVE_STEADY_CORRECTION_H
#define NEKTAR_WJ_NATIVE_STEADY_CORRECTION_H

#include <IncNavierStokesSolver/EquationSystems/WJNativeSteadyKrylov.h>
#include <Eigen/SparseLU>

namespace Nektar::WJNative
{
struct SteadySparseWork
{
    int colors=0, nonzeros=0, factorNonzeros=0, residualEvaluations=0, rejectedSamples=0;
    double baseNorm=0., assemblySeconds=0., factorSeconds=0., maximumUnassignedAction=0.;
    double inverseMass=0.,stiffness=0.;
    double fdRelativeStep=0.;
    std::string lastRejection,baseSha256;
};

/** A row-disjoint, centered finite-difference inverse on the free unknowns.
 * Supports describe the nonlinear rows of each column. The caller supplies
 * the exact linear entries on disjoint row/column pairs, such as the mixed
 * pressure and constraint action. The built matrix is Dr J Dx. Matrix-free
 * actions and checked linear residuals remain the responsibility of GMRES.
 */
class SteadySparseInverse
{
public:
    using SparseMatrix=Eigen::SparseMatrix<double>;
    using Action=std::function<SteadyVector(const SteadyVector&)>;

    SteadySparseInverse(std::vector<std::vector<int>> rowsByColumn,
        SparseMatrix exact, SteadyVector unknownScale, SteadyVector rowScale);
    void Prepare(const SteadyVector &base,const Action &residual,const Action &correction,
        const SteadyKrylovOptions &options,const std::function<void()> &deadline,SteadySparseWork &work);
    SteadyVector OriginalAction(const SteadyVector &direction)const;
    SteadyColoredColumnProbe ColumnProbe(int column)const;
    SteadyVector Solve(const SteadyVector &rhs)const;

private:
    std::vector<std::vector<int>> m_support,m_colors;
    SparseMatrix m_exact,m_original;
    SteadyVector m_unknownScale,m_rowScale;
    Eigen::SparseLU<SparseMatrix,Eigen::COLAMDOrdering<int>> m_factor;
};
}
#endif
