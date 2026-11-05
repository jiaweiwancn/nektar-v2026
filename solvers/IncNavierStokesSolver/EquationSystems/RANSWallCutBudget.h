///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_RANS_WALL_CUT_BUDGET_H
#define NEKTAR_RANS_WALL_CUT_BUDGET_H
namespace Nektar::RANSWallAlgebra
{
// A cut admission changes the working set. Even the last admitted cut must
// receive a subsequent reconstruction/oracle check, without another admission.
struct CutBudgetR12
{
    int limit;
    constexpr int Scans() const { return limit+1; }
    constexpr bool CanAdmit(int admitted) const { return admitted<limit; }
};
}
#endif
