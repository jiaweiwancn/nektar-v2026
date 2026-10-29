///////////////////////////////////////////////////////////////////////////////
// File: HITPlaneEvaluator.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Turbulence/HITPlaneEvaluator.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <new>
#include <stdexcept>
#include <utility>

namespace Nektar::SolverUtils
{
namespace
{
using Complex = std::complex<NekDouble>;

// Match the reference generator's rounded period and long-double reduction.
const NekDouble twoPi = 2.0 * std::acos(-1.0);

Complex Phase(long double value)
{
    const auto reduced = std::remainder(value,
                                       static_cast<long double>(twoPi));
    if (!std::isfinite(reduced))
    {
        throw std::overflow_error("HIT plane evaluator: phase overflow");
    }
    return {std::cos(static_cast<NekDouble>(reduced)),
            std::sin(static_cast<NekDouble>(reduced))};
}

bool AddPayload(std::size_t count, std::size_t size, std::size_t cap,
                std::size_t &total)
{
    if (total > cap || (size && count > (cap - total) / size))
    {
        return false;
    }
    total += count * size;
    return true;
}

NekDouble Velocity(long double value)
{
    if (!std::isfinite(value) ||
        std::abs(value) > std::numeric_limits<NekDouble>::max())
    {
        throw std::overflow_error("HIT plane evaluator: sample velocity overflow");
    }
    return static_cast<NekDouble>(value);
}
} // namespace

HITPlaneEvaluator::HITPlaneEvaluator(
    const HITGenerator &generator, const std::vector<HITVector> &points,
    std::size_t maxCacheBytes)
    : m_generator(generator), m_points(points)
{
    for (const auto &point : m_points)
    {
        for (const auto coordinate : point)
        {
            if (!std::isfinite(coordinate))
            {
                throw std::invalid_argument(
                    "HIT plane evaluator: coordinates must be finite");
            }
        }
    }
    if (m_points.empty())
    {
        m_grouped = true;
        m_method = "EmptyPlane";
        return;
    }
    for (const auto &point : m_points)
    {
        if (point[0] != m_points.front()[0])
        {
            m_method = "ReferenceNonplanar";
            return;
        }
    }
    const auto &modes = generator.GetModes();
    if (modes.empty())
    {
        m_grouped = true;
        m_method = "MeanOnly";
        return;
    }
    for (const auto &mode : modes)
    {
        for (int direction = 0; direction < 3; ++direction)
        {
            m_maxAbsWave[direction] = std::max(
                m_maxAbsWave[direction], std::abs(mode.wavevector[direction]));
        }
    }
    const auto &config = generator.GetConfig();
    if (!SafeAtTime(config.timeOrigin) ||
        !SafeAtTime(config.timeOrigin + config.recordDuration))
    {
        m_method = "ReferenceLargePhase";
        return;
    }

    std::size_t termBytes = 0;
    if (!AddPayload(modes.size(), sizeof(Term), maxCacheBytes, termBytes))
    {
        m_method = "ReferenceCacheLimit";
        return;
    }

    try
    {
        std::map<int, std::size_t> xIndex;
        std::map<std::pair<int, int>, std::size_t> transverseIndex;
        std::vector<Term> terms;
        std::vector<NekDouble> xWaves;
        std::vector<std::array<NekDouble, 2>> transverseWaves;
        terms.reserve(modes.size());
        auto bin = [&](int y, int z, NekDouble ky, NekDouble kz)
        {
            const auto entry = transverseIndex.emplace(
                std::make_pair(y, z), transverseWaves.size());
            if (entry.second)
            {
                transverseWaves.push_back({ky, kz});
            }
            return entry.first->second;
        };

        for (const auto &mode : modes)
        {
            const auto x = xIndex.emplace(mode.index[0], xWaves.size());
            if (x.second)
            {
                xWaves.push_back(mode.wavevector[0]);
            }
            const auto positive = bin(mode.index[1], mode.index[2],
                                      mode.wavevector[1], mode.wavevector[2]);
            const auto negative = bin(-mode.index[1], -mode.index[2],
                                      -mode.wavevector[1], -mode.wavevector[2]);
            terms.push_back({&mode, x.first->second, positive, negative});
        }

        std::size_t bytes = 0;
        const auto bins = transverseWaves.size();
        // Account for both persistent cache and Evaluate's temporary vectors.
        if (!AddPayload(terms.capacity(), sizeof(Term), maxCacheBytes, bytes) ||
            !AddPayload(xWaves.capacity(), sizeof(NekDouble),
                        maxCacheBytes, bytes) ||
            !AddPayload(xWaves.size(), sizeof(Complex), maxCacheBytes, bytes) ||
            !AddPayload(transverseWaves.capacity(),
                        sizeof(std::array<NekDouble, 2>), maxCacheBytes, bytes) ||
            !AddPayload(bins, sizeof(HITComplexVector), maxCacheBytes, bytes) ||
            (bins && m_points.size() >
                         std::numeric_limits<std::size_t>::max() / bins) ||
            !AddPayload(m_points.size() * bins, sizeof(Complex),
                        maxCacheBytes, bytes))
        {
            m_method = "ReferenceCacheLimit";
            return;
        }

        std::vector<Complex> spatialPhases(m_points.size() * bins);
        // The standard permits capacity to exceed size even for construction
        // with a count; retain the cap if an implementation does so.
        if (spatialPhases.capacity() > spatialPhases.size() &&
            !AddPayload(spatialPhases.capacity() - spatialPhases.size(),
                        sizeof(Complex), maxCacheBytes, bytes))
        {
            m_method = "ReferenceCacheLimit";
            return;
        }
        for (std::size_t point = 0; point < m_points.size(); ++point)
        {
            const long double dy =
                static_cast<long double>(m_points[point][1]) - config.origin[1];
            const long double dz =
                static_cast<long double>(m_points[point][2]) - config.origin[2];
            for (std::size_t index = 0; index < bins; ++index)
            {
                spatialPhases[point * bins + index] = Phase(
                    transverseWaves[index][0] * dy +
                    transverseWaves[index][1] * dz);
            }
        }
        m_terms = std::move(terms);
        m_xWaves = std::move(xWaves);
        m_transverseWaves = std::move(transverseWaves);
        m_spatialPhases = std::move(spatialPhases);
        m_cacheBytes = bytes;
        m_grouped = true;
        m_method = "GroupedTransverse";
    }
    catch (const std::bad_alloc &)
    {
        m_method = "ReferenceAllocation";
    }
    catch (const std::length_error &)
    {
        m_method = "ReferenceCacheLimit";
    }
}

bool HITPlaneEvaluator::SafeAtTime(NekDouble time) const
{
    // A conservative absolute-phase bound protects separated argument reduction
    // on platforms where long double has either 53 or 64 mantissa bits. Large
    // component phases are rejected even if their sum happens to cancel.
    const long double limit = std::min(
        1e6L, 1e-14L / std::numeric_limits<long double>::epsilon());
    const auto &config = m_generator.GetConfig();
    for (const auto &point : m_points)
    {
        long double dx = static_cast<long double>(point[0]) - config.origin[0];
        dx -= static_cast<long double>(config.meanSpeed) *
              (static_cast<long double>(time) - config.timeOrigin);
        const long double dy =
            static_cast<long double>(point[1]) - config.origin[1];
        const long double dz =
            static_cast<long double>(point[2]) - config.origin[2];
        const long double magnitude = m_maxAbsWave[0] * std::abs(dx) +
                                      m_maxAbsWave[1] * std::abs(dy) +
                                      m_maxAbsWave[2] * std::abs(dz);
        if (!std::isfinite(magnitude) || magnitude > limit)
        {
            return false;
        }
    }
    return true;
}

void HITPlaneEvaluator::Reference(
    NekDouble time, std::vector<HITVector> &velocities) const
{
    velocities.resize(m_points.size());
    for (std::size_t point = 0; point < m_points.size(); ++point)
    {
        velocities[point] = m_generator.Evaluate(m_points[point], time).velocity;
    }
}

void HITPlaneEvaluator::Evaluate(
    NekDouble time, std::vector<HITVector> &velocities) const
{
    m_generator.ValidateTime(time);
    if (!m_grouped || !SafeAtTime(time))
    {
        Reference(time, velocities);
        return;
    }
    const auto &config = m_generator.GetConfig();
    velocities.assign(m_points.size(), HITVector{config.meanSpeed, 0.0, 0.0});
    if (m_points.empty() || m_terms.empty())
    {
        return;
    }

    long double dx = static_cast<long double>(m_points.front()[0]) -
                     config.origin[0];
    dx -= static_cast<long double>(config.meanSpeed) *
          (static_cast<long double>(time) - config.timeOrigin);
    std::vector<Complex> xPhases(m_xWaves.size());
    for (std::size_t index = 0; index < m_xWaves.size(); ++index)
    {
        xPhases[index] = Phase(m_xWaves[index] * dx);
    }
    std::vector<HITComplexVector> amplitude(m_transverseWaves.size());
    for (const auto &term : m_terms)
    {
        for (int component = 0; component < 3; ++component)
        {
            const auto value = term.mode->coefficient[component] *
                               xPhases[term.xPhase];
            if (term.positiveBin == term.negativeBin)
            {
                amplitude[term.positiveBin][component] +=
                    Complex(2.0 * value.real(), 0.0);
            }
            else
            {
                amplitude[term.positiveBin][component] += value;
                amplitude[term.negativeBin][component] += std::conj(value);
            }
        }
    }
    const auto bins = amplitude.size();
    for (std::size_t point = 0; point < m_points.size(); ++point)
    {
        std::array<long double, 3> sum{config.meanSpeed, 0.0L, 0.0L};
        for (std::size_t index = 0; index < bins; ++index)
        {
            const auto spatial = m_spatialPhases[point * bins + index];
            for (int component = 0; component < 3; ++component)
            {
                sum[component] += (amplitude[index][component] * spatial).real();
            }
        }
        for (int component = 0; component < 3; ++component)
        {
            velocities[point][component] = Velocity(sum[component]);
        }
    }
}

std::size_t HITPlaneEvaluator::CacheBytes() const
{
    return m_cacheBytes;
}

bool HITPlaneEvaluator::IsGrouped() const
{
    return m_grouped;
}

const std::string &HITPlaneEvaluator::Method() const
{
    return m_method;
}
} // namespace Nektar::SolverUtils
