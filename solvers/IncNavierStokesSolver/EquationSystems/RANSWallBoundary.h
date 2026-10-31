///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_RANS_WALL_BOUNDARY_H
#define NEKTAR_RANS_WALL_BOUNDARY_H
#include <MultiRegions/ExpList.h>

namespace Nektar
{
class RANSWallBoundary
{
public:
    using Fields = Array<OneD, Array<OneD, NekDouble>>;
    using Gradients = Array<OneD, Fields>;
    virtual ~RANSWallBoundary() = default;
    virtual void CheckInitialState(const Fields &) const = 0;
    virtual std::string RestartSignature() const = 0;
    virtual void Update(const Fields &, const Fields &, const Gradients &,
                        const Fields &, const Fields &, const Fields &, NekDouble) = 0;
    virtual void AddRobinAuditResidual(int, const Fields &, Array<OneD, NekDouble> &) const = 0;
    virtual void WritePointDiagnostics(NekDouble, int) = 0;
    virtual void WriteDiagnostics(NekDouble, const Fields &, const Fields &,
                                  const Fields &, const Array<OneD, NekDouble> &) = 0;
    virtual NekDouble ConstraintResidual() const { return 0; }
    virtual Array<OneD,NekDouble> IterationState() const { return {}; }
    virtual void SetIterationState(const Array<OneD,NekDouble> &) {}
    virtual Array<OneD,NekDouble> IterationMetric() const { return {}; }
    virtual NekDouble IterationReactionScale() const { return 1; }
    virtual void SetTrialMode(bool) {}
};
}
#endif
