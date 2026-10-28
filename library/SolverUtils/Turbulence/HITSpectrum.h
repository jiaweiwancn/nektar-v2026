///////////////////////////////////////////////////////////////////////////////
// File: HITSpectrum.h
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_HITSPECTRUM_H
#define NEKTAR_SOLVERUTILS_HITSPECTRUM_H

#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <SolverUtils/SolverUtilsDeclspec.h>
#include <istream>
#include <string>
#include <vector>

namespace Nektar::SolverUtils
{
enum class HITInterpolation
{
    Linear,
    LogLinear
};

enum class HITLowWavenumber
{
    Zero,
    K4
};

enum class HITHighWavenumber
{
    Zero,
    Gaussian
};

/// Nonnegative shell spectrum: integral E(k) dk is turbulent kinetic energy.
/// Inputs and outputs use SI units and angular wavenumber. Invalid inputs throw
/// std::invalid_argument; arithmetic range failures throw std::overflow_error.
class SOLVER_UTILS_EXPORT HITSpectrum
{
public:
    static HITSpectrum VonKarman(NekDouble length, NekDouble amplitude);
    static HITSpectrum Tabulated(
        const std::vector<NekDouble> &wavenumbers,
        const std::vector<NekDouble> &energies,
        HITInterpolation interpolation = HITInterpolation::LogLinear,
        HITLowWavenumber low = HITLowWavenumber::K4,
        HITHighWavenumber high = HITHighWavenumber::Gaussian);

    /// First nonblank line must declare Shell3D, rad/m, m3/s2 exactly as in the
    /// version-1 contract. No numerical tail rows are inserted into the table.
    static HITSpectrum ReadTable(
        std::istream &input,
        HITInterpolation interpolation = HITInterpolation::LogLinear,
        HITLowWavenumber low = HITLowWavenumber::K4,
        HITHighWavenumber high = HITHighWavenumber::Gaussian);
    static HITSpectrum FromFile(
        const std::string &filename,
        HITInterpolation interpolation = HITInterpolation::LogLinear,
        HITLowWavenumber low = HITLowWavenumber::K4,
        HITHighWavenumber high = HITHighWavenumber::Gaussian);

    NekDouble Evaluate(NekDouble k) const;
    bool IsIdenticallyZero() const;
    const std::vector<NekDouble> &GetWavenumbers() const;
    const std::vector<NekDouble> &GetEnergies() const;

private:
    HITSpectrum() = default;
    bool m_tabulated = false;
    bool m_identicallyZero = false;
    NekDouble m_length = 1.0;
    NekDouble m_amplitude = 1.0;
    std::vector<NekDouble> m_wavenumbers;
    std::vector<NekDouble> m_energies;
    HITInterpolation m_interpolation = HITInterpolation::LogLinear;
    HITLowWavenumber m_low = HITLowWavenumber::K4;
    HITHighWavenumber m_high = HITHighWavenumber::Gaussian;
    NekDouble m_gaussianLogBeta = 0.0;
};
} // namespace Nektar::SolverUtils

#endif
