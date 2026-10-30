///////////////////////////////////////////////////////////////////////////////
// File: RANSHistory.h
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_RANS_HISTORY_H
#define NEKTAR_RANS_HISTORY_H
#include <LibUtilities/TimeIntegration/IMEXTimeIntegrationSchemes.h>
#include <MultiRegions/ExpList.h>

namespace Nektar
{
/** Native IMEX tableaux with complete quadrature-history restoration. The
 * underlying algorithms, stage coefficients and operators are unchanged. */
class RANSIMEX : public LibUtilities::IMEXTimeIntegrationScheme
{
public:
    explicit RANSIMEX(size_t order) : IMEXTimeIntegrationScheme("", order, {}) {}
    void Restore(const LibUtilities::TripleArray &values,
                 const Array<OneD, NekDouble> &times, size_t steps);
private:
    void v_InitializeScheme(NekDouble dt, LibUtilities::ConstDoubleArray &initial,
        NekDouble time, const LibUtilities::TimeIntegrationSchemeOperators &op) override;
    LibUtilities::ConstDoubleArray &v_TimeIntegrate(size_t step, NekDouble dt) override;
    LibUtilities::TripleArray m_restored;
    Array<OneD, NekDouble> m_times;
    size_t m_offset = 0;
};

struct RANSHistory
{
    using Fields = Array<OneD, MultiRegions::ExpListSharedPtr>;
    static std::string Write(const std::string &directory, const Fields &fields,
                      const RANSIMEX &scheme, const std::string &signature,
                      size_t steps, NekDouble time, NekDouble dt);
    static size_t Read(const std::string &directory, const Fields &fields,
                       const LibUtilities::DoubleArray &importedCoefficients,
                       RANSIMEX &scheme, const std::string &signature,
                       NekDouble time, NekDouble dt);
};
}
#endif
