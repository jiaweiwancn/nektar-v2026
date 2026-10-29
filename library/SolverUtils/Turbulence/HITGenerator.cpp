///////////////////////////////////////////////////////////////////////////////
// File: HITGenerator.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Turbulence/HITGenerator.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace Nektar::SolverUtils
{
namespace
{
const NekDouble pi = std::acos(-1.0);
const NekDouble twoPi = 2.0 * pi;
using WideTensor = std::array<std::array<long double, 3>, 3>;

void Require(bool condition, const std::string &message)
{
    if (!condition)
    {
        throw std::invalid_argument("HIT generator: " + message);
    }
}

NekDouble Narrow(long double value, const std::string &quantity)
{
    if (!std::isfinite(value) ||
        std::abs(value) > std::numeric_limits<NekDouble>::max())
    {
        throw std::overflow_error("HIT generator: unrepresentable " + quantity);
    }
    return static_cast<NekDouble>(value);
}

void ValidateCounts(const HITIndex &counts)
{
    // 6*(product-1)+5 must fit before the counter is evaluated. Check the
    // multiplication itself, rather than detecting wrap after multiplication.
    const auto limit = std::numeric_limits<std::uint64_t>::max() / 6;
    std::uint64_t product = 1;
    for (int count : counts)
    {
        Require(count >= 4 && count % 2 == 0,
                "grid counts must be even integers at least four");
        if (product > limit / static_cast<std::uint64_t>(count))
        {
            throw std::overflow_error("HIT generator: mode counter overflow");
        }
        product *= static_cast<std::uint64_t>(count);
    }
}

struct Lattice
{
    HITVector spacing;
    HITIndex limits;
    long double volume;
};

Lattice Validate(const HITGeneratorConfig &config)
{
    ValidateCounts(config.counts);
    Require(std::isfinite(config.meanSpeed) && config.meanSpeed > 0.0,
            "mean speed must be finite and positive");
    Require(std::isfinite(config.timeStep) && config.timeStep > 0.0,
            "time step must be finite and positive");
    Require(std::isfinite(config.timeOrigin), "time origin must be finite");
    Require(std::isfinite(config.recordDuration) && config.recordDuration > 0.0,
            "record duration must be finite and positive");
    Require(std::isfinite(config.kMax) && config.kMax > 0.0 &&
                std::isfinite(config.taperStart) && config.taperStart >= 0.0 &&
                config.taperStart < config.kMax,
            "require 0 <= taper start < positive cutoff");
    Require(config.normalization == HITNormalization::Preserve ||
                config.normalization == HITNormalization::TargetRMS,
            "unknown normalization");
    Require(std::isfinite(config.targetRMS) && config.targetRMS >= 0.0,
            "target RMS must be finite and nonnegative");
    const NekDouble end = Narrow(
        static_cast<long double>(config.timeOrigin) + config.recordDuration,
        "record end time");
    Require(end > config.timeOrigin,
            "record duration is not representable at the time origin");
    const long double period =
        static_cast<long double>(config.lengths[0]) / config.meanSpeed;
    Require(config.allowRepeat || config.recordDuration < period,
            "record duration must be shorter than the synthetic period");
    Require(static_cast<long double>(config.meanSpeed) * config.kMax *
                    config.timeStep < pi,
            "temporal Nyquist bound exceeded");

    Lattice lattice;
    lattice.volume = 1.0L;
    for (std::size_t d = 0; d < 3; ++d)
    {
        Require(std::isfinite(config.origin[d]), "spatial origin must be finite");
        Require(std::isfinite(config.lengths[d]) && config.lengths[d] > 0.0,
                "box periods must be finite and positive");
        lattice.spacing[d] = Narrow(
            static_cast<long double>(twoPi) / config.lengths[d],
            "wavenumber spacing");
        Require(lattice.spacing[d] > 0.0, "wavenumber spacing underflow");
        const long double nyquist =
            static_cast<long double>(lattice.spacing[d]) * (config.counts[d] / 2);
        Require(config.kMax < nyquist,
                "radial cutoff must be below every Nyquist plane");
        lattice.volume *= lattice.spacing[d];
        const auto radialLimit = static_cast<int>(
            std::floor(config.kMax / lattice.spacing[d]));
        lattice.limits[d] =
            std::min(config.counts[d] / 2 - 1, radialLimit);
    }
    Require(std::isfinite(lattice.volume) && lattice.volume > 0.0L,
            "spectral cell volume is not representable");
    return lattice;
}

bool MeanOnly(const HITGeneratorConfig &config, const HITSpectrum &spectrum)
{
    return (config.normalization == HITNormalization::TargetRMS &&
            config.targetRMS == 0.0) ||
           (config.normalization == HITNormalization::Preserve &&
            spectrum.IsIdenticallyZero());
}

template <typename Callback>
void Visit(const HITGeneratorConfig &config, const HITSpectrum &spectrum,
           const Lattice &lattice, Callback callback)
{
    // The first nonzero integer component is positive. Bounds avoid visiting
    // the negative conjugate half, Nyquist planes, and the zero vector.
    for (int mx = 0; mx <= lattice.limits[0]; ++mx)
    {
        for (int my = mx == 0 ? 0 : -lattice.limits[1];
             my <= lattice.limits[1]; ++my)
        {
            const int firstZ = mx == 0 && my == 0 ? 1 : -lattice.limits[2];
            for (int mz = firstZ; mz <= lattice.limits[2]; ++mz)
            {
                const HITIndex index{mx, my, mz};
                HITVector wave{}, unit{};
                for (std::size_t d = 0; d < 3; ++d)
                {
                    wave[d] = index[d] * lattice.spacing[d];
                }
                const NekDouble k = std::hypot(wave[0], wave[1], wave[2]);
                if (k >= config.kMax)
                {
                    continue;
                }
                NekDouble energy = spectrum.Evaluate(k);
                if (k > config.taperStart)
                {
                    // cos^2 of the half angle avoids cancellation close to
                    // the cutoff; this is the contract's energy taper W(k).
                    const NekDouble taper = std::cos(
                        0.5 * pi * (k - config.taperStart) /
                        (config.kMax - config.taperStart));
                    energy *= taper * taper;
                }
                if (energy == 0.0)
                {
                    continue;
                }
                for (std::size_t d = 0; d < 3; ++d)
                {
                    unit[d] = wave[d] / k;
                }
                // Phi = weight * (I - unit*unit^T). Avoid forming k^4.
                const long double weight =
                    lattice.volume * energy /
                    (4.0L * pi * static_cast<long double>(k) * k);
                if (!std::isfinite(weight) || weight < 0.0L)
                {
                    throw std::overflow_error("HIT generator: spectral weight overflow");
                }
                callback(index, wave, unit, k, weight);
            }
        }
    }
}

struct Summary
{
    std::size_t count = 0;
    std::size_t bytes = 0;
    long double energy = 0.0L;
    long double k2Energy = 0.0L;
    WideTensor covariance{};
};

Summary Preflight(const HITGeneratorConfig &config, const HITSpectrum &spectrum,
                  const Lattice &lattice)
{
    Summary result;
    if (MeanOnly(config, spectrum))
    {
        return result;
    }
    Visit(config, spectrum, lattice,
          [&](const HITIndex &, const HITVector &, const HITVector &unit,
              NekDouble k, long double weight)
          {
              if (result.count ==
                  std::numeric_limits<std::size_t>::max() / sizeof(HITFourierMode))
              {
                  throw std::overflow_error("HIT generator: mode allocation overflow");
              }
              ++result.count;
              result.bytes = result.count * sizeof(HITFourierMode);
              if (result.bytes > config.maxModeBytes)
              {
                  throw std::length_error("HIT generator: mode memory limit exceeded");
              }
              // Full-lattice K is twice the scalar weight per independent
              // member. Full covariance includes both conjugate partners.
              result.energy += 2.0L * weight;
              result.k2Energy += 2.0L * weight * k * k;
              for (std::size_t i = 0; i < 3; ++i)
              {
                  for (std::size_t j = 0; j < 3; ++j)
                  {
                      result.covariance[i][j] +=
                          2.0L * weight * ((i == j ? 1.0L : 0.0L) -
                                          static_cast<long double>(unit[i]) * unit[j]);
                  }
              }
          });
    if (config.normalization == HITNormalization::TargetRMS)
    {
        Require(result.energy > 0.0L,
                "positive target RMS has no represented spectral energy");
    }
    if (!std::isfinite(result.energy) || !std::isfinite(result.k2Energy))
    {
        throw std::overflow_error("HIT generator: spectral moments overflow");
    }
    return result;
}
} // namespace

HITGenerator::HITGenerator(const HITGeneratorConfig &config,
                         const HITSpectrum &spectrum)
    : m_config(config)
{
    const auto lattice = Validate(config);
    const auto summary = Preflight(config, spectrum, lattice);
    if (MeanOnly(config, spectrum))
    {
        m_diagnostics.energyScale =
            config.normalization == HITNormalization::TargetRMS ? 0.0 : 1.0;
        return;
    }
    const long double requestedScale =
        config.normalization == HITNormalization::TargetRMS
            ? 1.5L * config.targetRMS * config.targetRMS / summary.energy
            : 1.0L;
    m_diagnostics.energyScale = Narrow(requestedScale, "energy normalization");
    if (m_diagnostics.energyScale == 0.0 ||
        std::abs(static_cast<long double>(m_diagnostics.energyScale) /
                     requestedScale - 1.0L) >
            64.0L * std::numeric_limits<NekDouble>::epsilon())
    {
        throw std::overflow_error(
            "HIT generator: energy normalization underflows numerical precision");
    }
    // The stored normalization is the one actually used for coefficients and
    // diagnostics, including on subsequent reproducible reconstructions.
    const long double scale = m_diagnostics.energyScale;
    m_diagnostics.independentModes = summary.count;
    m_diagnostics.modeBytes = summary.bytes;
    m_diagnostics.expectedKineticEnergy =
        Narrow(summary.energy * scale, "expected energy");
    m_expectedK2Energy = Narrow(summary.k2Energy * scale, "expected k-squared energy");
    for (std::size_t i = 0; i < 3; ++i)
    {
        for (std::size_t j = 0; j < 3; ++j)
        {
            m_diagnostics.expectedCovariance[i][j] =
                Narrow(summary.covariance[i][j] * scale, "expected covariance");
        }
    }

    m_modes.reserve(summary.count);
    WideTensor realizedCovariance{};
    long double realizedEnergy = 0.0L, realizedK2 = 0.0L;
    Visit(config, spectrum, lattice,
          [&](const HITIndex &index, const HITVector &wave,
              const HITVector &unit, NekDouble k, long double weight)
          {
              const auto xi = Gaussian(config.seed, ModeId(index, config.counts));
              const NekDouble factor = Narrow(std::sqrt(weight * scale), "mode amplitude");
              HITFourierMode mode{index, wave, {}};
              mode.coefficient[0] = factor * (xi[1] * unit[2] - xi[2] * unit[1]);
              mode.coefficient[1] = factor * (xi[2] * unit[0] - xi[0] * unit[2]);
              mode.coefficient[2] = factor * (xi[0] * unit[1] - xi[1] * unit[0]);
              long double modeEnergy = 0.0L;
              for (std::size_t i = 0; i < 3; ++i)
              {
                  const auto ai = mode.coefficient[i];
                  modeEnergy += static_cast<long double>(ai.real()) * ai.real() +
                                static_cast<long double>(ai.imag()) * ai.imag();
                  for (std::size_t j = 0; j < 3; ++j)
                  {
                      const auto aj = mode.coefficient[j];
                      realizedCovariance[i][j] +=
                          2.0L * (static_cast<long double>(ai.real()) * aj.real() +
                                  static_cast<long double>(ai.imag()) * aj.imag());
                  }
              }
              realizedEnergy += modeEnergy;
              realizedK2 += modeEnergy * k * k;
              m_modes.push_back(mode);
          });
    m_diagnostics.realizedKineticEnergy = Narrow(realizedEnergy, "realized energy");
    m_realizedK2Energy = Narrow(realizedK2, "realized k-squared energy");
    for (std::size_t i = 0; i < 3; ++i)
    {
        for (std::size_t j = 0; j < 3; ++j)
        {
            m_diagnostics.realizedCovariance[i][j] =
                Narrow(realizedCovariance[i][j], "realized covariance");
        }
    }
}

const HITGeneratorConfig &HITGenerator::GetConfig() const
{
    return m_config;
}

const std::vector<HITFourierMode> &HITGenerator::GetModes() const
{
    return m_modes;
}

std::size_t HITGenerator::EstimateModeBytes(const HITGeneratorConfig &config,
                                          const HITSpectrum &spectrum)
{
    const auto lattice = Validate(config);
    return Preflight(config, spectrum, lattice).bytes;
}

void HITGenerator::ValidateTime(NekDouble time) const
{
    Require(std::isfinite(time), "evaluation time must be finite");
    const NekDouble end = m_config.timeOrigin + m_config.recordDuration;
    const NekDouble tolerance = std::min(
        64.0 * std::numeric_limits<NekDouble>::epsilon() *
            std::max({1.0, std::abs(m_config.timeOrigin), std::abs(end)}),
        1e-10 * m_config.recordDuration);
    Require(time >= m_config.timeOrigin - tolerance && time <= end + tolerance,
            "evaluation time is outside the declared record");
}

HITSample HITGenerator::Evaluate(const HITVector &point, NekDouble time) const
{
    ValidateTime(time);
    std::array<long double, 3> displacement;
    for (std::size_t d = 0; d < 3; ++d)
    {
        Require(std::isfinite(point[d]), "evaluation coordinates must be finite");
        displacement[d] = static_cast<long double>(point[d]) - m_config.origin[d];
    }
    displacement[0] -= static_cast<long double>(m_config.meanSpeed) *
                       (static_cast<long double>(time) - m_config.timeOrigin);
    HITSample result;
    result.velocity[0] = m_config.meanSpeed;
    for (const auto &mode : m_modes)
    {
        long double phase = 0.0L;
        for (std::size_t d = 0; d < 3; ++d)
        {
            phase += mode.wavevector[d] * displacement[d];
        }
        phase = std::remainder(phase, static_cast<long double>(twoPi));
        if (!std::isfinite(phase))
        {
            throw std::overflow_error("HIT generator: phase overflow");
        }
        const std::complex<NekDouble> exponential(
            std::cos(static_cast<NekDouble>(phase)),
            std::sin(static_cast<NekDouble>(phase)));
        for (std::size_t i = 0; i < 3; ++i)
        {
            const auto value = mode.coefficient[i] * exponential;
            result.velocity[i] += 2.0 * value.real();
            for (std::size_t j = 0; j < 3; ++j)
            {
                result.gradient[i][j] -= 2.0 * mode.wavevector[j] * value.imag();
            }
            result.timeDerivative[i] +=
                2.0 * m_config.meanSpeed * mode.wavevector[0] * value.imag();
        }
    }
    for (std::size_t i = 0; i < 3; ++i)
    {
        Narrow(result.velocity[i], "sample velocity");
        Narrow(result.timeDerivative[i], "sample time derivative");
        for (std::size_t j = 0; j < 3; ++j)
        {
            Narrow(result.gradient[i][j], "sample gradient");
        }
    }
    return result;
}

HITDiagnostics HITGenerator::Diagnostics(NekDouble viscosity) const
{
    Require(std::isfinite(viscosity) && viscosity >= 0.0,
            "diagnostic viscosity must be finite and nonnegative");
    auto result = m_diagnostics;
    result.expectedDissipation = Narrow(
        2.0L * viscosity * m_expectedK2Energy, "expected dissipation");
    result.realizedDissipation = Narrow(
        2.0L * viscosity * m_realizedK2Energy, "realized dissipation");
    return result;
}

std::uint64_t HITGenerator::ParseSeed(const std::string &text)
{
    Require(!text.empty(), "seed must be unsigned decimal text");
    std::uint64_t value = 0;
    for (char digit : text)
    {
        Require(digit >= '0' && digit <= '9', "seed must be unsigned decimal text");
        const auto d = static_cast<std::uint64_t>(digit - '0');
        if (value > (std::numeric_limits<std::uint64_t>::max() - d) / 10)
        {
            throw std::overflow_error("HIT generator: seed exceeds uint64 range");
        }
        value = 10 * value + d;
    }
    return value;
}

std::uint64_t HITGenerator::Mix64(std::uint64_t value)
{
    auto z = value + UINT64_C(0x9e3779b97f4a7c15);
    z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31);
}

NekDouble HITGenerator::Uniform52(std::uint64_t word)
{
    return (static_cast<NekDouble>(word >> 12) + 0.5) * 0x1p-52;
}

std::uint64_t HITGenerator::ModeId(const HITIndex &index, const HITIndex &counts)
{
    ValidateCounts(counts);
    bool canonical = false;
    for (std::size_t d = 0; d < 3; ++d)
    {
        Require(index[d] >= -counts[d] / 2 + 1 && index[d] <= counts[d] / 2 - 1,
                "mode index reaches an excluded Nyquist plane");
        if (!canonical && index[d] != 0)
        {
            Require(index[d] > 0, "mode must be the canonical conjugate member");
            canonical = true;
        }
    }
    Require(canonical, "the zero mode has no random coefficient");
    std::uint64_t id = 0;
    for (std::size_t d = 0; d < 3; ++d)
    {
        id = id * static_cast<std::uint64_t>(counts[d]) +
             static_cast<std::uint64_t>(index[d] + counts[d] / 2);
    }
    return id;
}

std::array<std::uint64_t, 6> HITGenerator::CounterWords(std::uint64_t seed,
                                                      std::uint64_t modeId)
{
    if (modeId > (std::numeric_limits<std::uint64_t>::max() - 5) / 6)
    {
        throw std::overflow_error("HIT generator: random counter overflow");
    }
    std::array<std::uint64_t, 6> result;
    for (std::uint64_t q = 0; q < result.size(); ++q)
    {
        result[q] = Mix64(seed + 6 * modeId + q);
    }
    return result;
}

HITComplexVector HITGenerator::Gaussian(std::uint64_t seed, std::uint64_t modeId)
{
    const auto words = CounterWords(seed, modeId);
    HITComplexVector result;
    for (std::size_t j = 0; j < 3; ++j)
    {
        const NekDouble radius = std::sqrt(-std::log(Uniform52(words[2 * j])));
        const NekDouble angle = twoPi * Uniform52(words[2 * j + 1]);
        // Radius already includes the complex 1/sqrt(2) normalization.
        result[j] = {radius * std::cos(angle), radius * std::sin(angle)};
    }
    return result;
}

const char *HITGenerator::Algorithm()
{
    return "SplitMix64BoxMuller52V1";
}
} // namespace Nektar::SolverUtils
