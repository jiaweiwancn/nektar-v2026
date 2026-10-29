///////////////////////////////////////////////////////////////////////////////
// File: HITGenerator.h
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_HITGENERATOR_H
#define NEKTAR_SOLVERUTILS_HITGENERATOR_H

#include <SolverUtils/Turbulence/HITSpectrum.h>
#include <array>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace Nektar::SolverUtils
{
using HITVector = std::array<NekDouble, 3>;
using HITIndex = std::array<int, 3>;
using HITComplexVector = std::array<std::complex<NekDouble>, 3>;
using HITTensor = std::array<HITVector, 3>;

enum class HITNormalization
{
    Preserve,
    TargetRMS
};

struct HITGeneratorConfig
{
    HITVector lengths{6.2831853071795864769, 6.2831853071795864769,
                      6.2831853071795864769};
    HITIndex counts{8, 8, 8};
    HITVector origin{0.0, 0.0, 0.0};
    std::uint64_t seed = 0;
    NekDouble meanSpeed = 1.0;
    NekDouble timeOrigin = 0.0;
    NekDouble recordDuration = 0.5;
    NekDouble timeStep = 0.001;
    NekDouble kMax = 3.5;
    NekDouble taperStart = 3.0;
    bool allowRepeat = false;
    HITNormalization normalization = HITNormalization::Preserve;
    NekDouble targetRMS = 0.0;
    /// A resource guard for native callers, independent of the realization.
    std::size_t maxModeBytes = std::numeric_limits<std::size_t>::max();
};

struct HITFourierMode
{
    HITIndex index;
    HITVector wavevector;
    HITComplexVector coefficient;
};

struct HITSample
{
    HITVector velocity{};
    /// gradient[i][j] = partial velocity component i / partial coordinate j.
    HITTensor gradient{};
    HITVector timeDerivative{};
};

struct HITDiagnostics
{
    std::size_t independentModes = 0;
    std::size_t modeBytes = 0;
    NekDouble energyScale = 1.0;
    NekDouble expectedKineticEnergy = 0.0;
    NekDouble realizedKineticEnergy = 0.0;
    NekDouble expectedDissipation = 0.0;
    NekDouble realizedDissipation = 0.0;
    HITTensor expectedCovariance{};
    HITTensor realizedCovariance{};
};

/// Reference three-dimensional Fourier synthesis. Stores canonical half-lattice
/// coefficients; Evaluate includes both conjugate partners and the mean flow.
/// This core has no solver, XML, MPI, boundary projection or FFT dependency.
class SOLVER_UTILS_EXPORT HITGenerator
{
public:
    HITGenerator(const HITGeneratorConfig &config, const HITSpectrum &spectrum);
    const HITGeneratorConfig &GetConfig() const;
    const std::vector<HITFourierMode> &GetModes() const;
    HITSample Evaluate(const HITVector &point, NekDouble time) const;
    HITDiagnostics Diagnostics(NekDouble viscosity = 0.0) const;

    /// Performs validation and spectral preflight without allocating modes or
    /// drawing random numbers; callers can report memory before construction.
    static std::size_t EstimateModeBytes(const HITGeneratorConfig &config,
                                         const HITSpectrum &spectrum);
    static std::uint64_t ParseSeed(const std::string &text);
    static std::uint64_t Mix64(std::uint64_t value);
    static NekDouble Uniform52(std::uint64_t word);
    static std::uint64_t ModeId(const HITIndex &index, const HITIndex &counts);
    static std::array<std::uint64_t, 6> CounterWords(std::uint64_t seed,
                                                   std::uint64_t modeId);
    static HITComplexVector Gaussian(std::uint64_t seed, std::uint64_t modeId);
    static const char *Algorithm();

private:
    friend class HITPlaneEvaluator;
    void ValidateTime(NekDouble time) const;
    HITGeneratorConfig m_config;
    std::vector<HITFourierMode> m_modes;
    HITDiagnostics m_diagnostics;
    NekDouble m_expectedK2Energy = 0.0;
    NekDouble m_realizedK2Energy = 0.0;
};
} // namespace Nektar::SolverUtils

#endif
