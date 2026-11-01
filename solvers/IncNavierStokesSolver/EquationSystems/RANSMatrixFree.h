///////////////////////////////////////////////////////////////////////////////
// File: RANSMatrixFree.h
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_RANS_MATRIX_FREE_H
#define NEKTAR_RANS_MATRIX_FREE_H
#include <MultiRegions/ContField.h>
#include <functional>

namespace Nektar
{
/** Current-coefficient native Helmholtz operator with a cached SPD
 * preconditioner. Only the preconditioner retains old coefficients. */
class RANSMatrixFree
{
public:
    using Fields = Array<OneD, Array<OneD, NekDouble>>;
    using RobinOperator = std::function<void(int, const Fields &, Array<OneD, NekDouble> &)>;
    RANSMatrixFree(const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                   int dimension, NekDouble tolerance, RobinOperator robin,
                   NekDouble gradDiv = 0);
    ~RANSMatrixFree();
    int Solve(int field, const Array<OneD, NekDouble> &rhs,
              const Array<OneD, NekDouble> &diffusion,
              const Array<OneD, NekDouble> &mass, NekDouble alpha,
              Array<OneD, NekDouble> &output);
    int SolveVelocity(const Fields &rhs, const Fields &diffusion, NekDouble alpha,
                      Fields &output);
private:
    struct Preconditioner
    {
        int field;
        NekDouble alpha;
        Array<OneD, NekDouble> diffusion, directional, mass;
        std::unique_ptr<MultiRegions::GlobalLinSysKey> key;
    };
    Array<OneD, MultiRegions::ExpListSharedPtr> m_fields;
    int m_dim;
    NekDouble m_tolerance, m_gradDiv;
    RobinOperator m_robin;
    std::vector<std::unique_ptr<Preconditioner>> m_preconditioners;
};
}
#endif
