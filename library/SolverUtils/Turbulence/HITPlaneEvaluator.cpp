///////////////////////////////////////////////////////////////////////////////
// File: HITPlaneEvaluator.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Turbulence/HITPlaneEvaluator.h>

#include <algorithm>
#include <cmath>
#include <limits>
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

bool Product(std::size_t a, std::size_t b, std::size_t &result)
{
    if (b && a > std::numeric_limits<std::size_t>::max() / b)
    {
        return false;
    }
    result = a * b;
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
    NekDouble maxCoefficient = 0.0;
    std::array<std::size_t, 3> maxIndex{};
    for (const auto &mode : modes)
    {
        for (int direction = 0; direction < 3; ++direction)
        {
            m_maxAbsWave[direction] = std::max(
                m_maxAbsWave[direction], std::abs(mode.wavevector[direction]));
            maxCoefficient = std::max(
                maxCoefficient, std::abs(mode.coefficient[direction]));
            const auto index = static_cast<std::int64_t>(mode.index[direction]);
            maxIndex[direction] = std::max(maxIndex[direction],
                static_cast<std::size_t>(index < 0 ? -index : index));
        }
    }
    const auto &config = generator.GetConfig();
    const NekDouble maxWave = std::max(
        {m_maxAbsWave[0], m_maxAbsWave[1], m_maxAbsWave[2]});
    const long double logWave = maxWave > 0.0
        ? std::log(static_cast<long double>(maxWave)) : 0.0L;
    const long double logSpeed = std::max(
        0.0L, std::log(static_cast<long double>(config.meanSpeed)));
    const long double logLimit = std::log(
        static_cast<long double>(std::numeric_limits<NekDouble>::max()) / 16.0L);
    // Reference() also evaluates gradients and time derivatives. Bound their
    // paired sums without forming products or using possibly underflowed K.
    // Separately protect the reference's intermediate 2*U*k multiplications,
    // even when the subsequent coefficient is tiny or zero.
    const long double logPartial =
        std::log(2.0L) + std::max(0.0L, logWave) + logSpeed;
    long double logDerivative = 0.0L;
    if (maxWave > 0.0 && maxCoefficient > 0.0)
    {
        logDerivative = std::log(2.0L) +
            std::log(static_cast<long double>(modes.size())) + logWave +
            std::log(static_cast<long double>(maxCoefficient)) + logSpeed;
    }
    if (!std::isfinite(logPartial) || !std::isfinite(logDerivative) ||
        logPartial > logLimit || logDerivative > logLimit)
    {
        m_method = "ReferenceDerivativeRange";
        return;
    }
    if (!SafeAtTime(config.timeOrigin) ||
        !SafeAtTime(config.timeOrigin + config.recordDuration))
    {
        m_method = "ReferenceLargePhase";
        return;
    }

    // Bound dense lattice lookup sizes before allocating them. Their storage
    // grows with the index extents, not with the number of retained 3D modes.
    // Use only occupied extents: session counts can greatly exceed the cutoff.
    const auto maxSize = std::numeric_limits<std::size_t>::max();
    if (maxIndex[0] == maxSize ||
        maxIndex[1] > (maxSize - 1) / 2 ||
        maxIndex[2] > (maxSize - 1) / 2)
    {
        m_method = "ReferenceCacheLimit";
        return;
    }
    const auto xSlots = maxIndex[0] + 1;
    const auto ySlots = 2 * maxIndex[1] + 1;
    const auto zSlots = 2 * maxIndex[2] + 1;
    std::size_t lookupBytes = 0;
    if (ySlots > maxSize / zSlots ||
        !AddPayload(xSlots, sizeof(std::size_t), maxCacheBytes, lookupBytes) ||
        !AddPayload(ySlots * zSlots, sizeof(std::size_t), maxCacheBytes,
                    lookupBytes))
    {
        m_method = "ReferenceCacheLimit";
        return;
    }

    try
    {
        const auto missing = std::numeric_limits<std::size_t>::max();
        std::vector<std::size_t> xIndex(xSlots, missing);
        std::vector<std::size_t> transverseIndex(ySlots * zSlots, missing);
        std::vector<NekDouble> xWaves;
        std::vector<std::array<NekDouble, 2>> transverseWaves;
        auto latticeIndex = [&](int y, int z)
        {
            const auto iy = static_cast<std::size_t>(
                static_cast<std::int64_t>(y) +
                static_cast<std::int64_t>(maxIndex[1]));
            const auto iz = static_cast<std::size_t>(
                static_cast<std::int64_t>(z) +
                static_cast<std::int64_t>(maxIndex[2]));
            return iy * zSlots + iz;
        };
        auto bin = [&](int y, int z, NekDouble ky, NekDouble kz)
        {
            auto &index = transverseIndex[latticeIndex(y, z)];
            if (index == missing)
            {
                index = transverseWaves.size();
                transverseWaves.push_back({ky, kz});
            }
        };

        for (const auto &mode : modes)
        {
            auto &x = xIndex[static_cast<std::size_t>(mode.index[0])];
            if (x == missing)
            {
                x = xWaves.size();
                xWaves.push_back(mode.wavevector[0]);
            }
            // Keep first-appearance bin order and conjugate insertion order
            // identical to the original per-mode grouped evaluator.
            bin(mode.index[1], mode.index[2],
                mode.wavevector[1], mode.wavevector[2]);
            bin(-mode.index[1], -mode.index[2],
                -mode.wavevector[1], -mode.wavevector[2]);
        }

        std::size_t bytes = 0;
        const auto bins = transverseWaves.size();
        // Account for both persistent cache and Evaluate's temporary vectors.
        if (!AddPayload(xIndex.capacity(), sizeof(std::size_t),
                        maxCacheBytes, bytes) ||
            !AddPayload(transverseIndex.capacity(), sizeof(std::size_t),
                        maxCacheBytes, bytes) ||
            !AddPayload(xWaves.capacity(), sizeof(NekDouble),
                        maxCacheBytes, bytes) ||
            !AddPayload(xWaves.size(), sizeof(Complex), maxCacheBytes, bytes) ||
            !AddPayload(transverseWaves.capacity(),
                        sizeof(std::array<NekDouble, 2>), maxCacheBytes, bytes) ||
            !AddPayload(bins, sizeof(HITComplexVector), maxCacheBytes, bytes))
        {
            m_method = "ReferenceCacheLimit";
            return;
        }

        const auto commonBytes = bytes;
        std::size_t phaseCount = 0;
        bool fullCache = Product(m_points.size(), bins, phaseCount) &&
            AddPayload(phaseCount, sizeof(Complex), maxCacheBytes, bytes);
        std::vector<Complex> spatialPhases;
        std::array<std::vector<Complex>, 2> axisPhases;
        std::vector<std::array<std::size_t, 2>> pointIndex;
        std::array<std::size_t, 2> uniqueCoordinates{};
        std::size_t firstAxis = 0;
        if (fullCache)
        {
            spatialPhases.resize(phaseCount);
            // The standard permits capacity to exceed size even for count
            // construction; retain the cap if an implementation does so.
            if (spatialPhases.capacity() > spatialPhases.size() &&
                !AddPayload(spatialPhases.capacity() - spatialPhases.size(),
                            sizeof(Complex), maxCacheBytes, bytes))
            {
                fullCache = false;
                std::vector<Complex>().swap(spatialPhases);
            }
        }
        if (fullCache)
        {
            for (std::size_t point = 0; point < m_points.size(); ++point)
            {
                const long double dy = static_cast<long double>(
                    m_points[point][1]) - config.origin[1];
                const long double dz = static_cast<long double>(
                    m_points[point][2]) - config.origin[2];
                for (std::size_t index = 0; index < bins; ++index)
                {
                    spatialPhases[point * bins + index] = Phase(
                        transverseWaves[index][0] * dy +
                        transverseWaves[index][1] * dz);
                }
            }
        }
        else
        {
            // Reuse one-dimensional phases when the full point-by-bin table
            // does not fit. Arbitrary point pairs and duplicates are mapped
            // to their unique y/z coordinates. The temporary coordinate lists
            // are discarded after constructing the bounded persistent cache.
            bytes = commonBytes;
            if (!AddPayload(m_points.size(), sizeof(pointIndex[0]),
                            maxCacheBytes, bytes))
            {
                m_method = "ReferenceCacheLimit";
                return;
            }
            pointIndex.resize(m_points.size());
            if (pointIndex.capacity() > pointIndex.size() &&
                !AddPayload(pointIndex.capacity() - pointIndex.size(),
                            sizeof(pointIndex[0]), maxCacheBytes, bytes))
            {
                m_method = "ReferenceCacheLimit";
                return;
            }
            std::array<std::vector<NekDouble>, 2> coordinates;
            std::array<std::vector<NekDouble>, 2> waves{
                std::vector<NekDouble>(ySlots), std::vector<NekDouble>(zSlots)};
            const std::array<std::size_t, 2> slots{ySlots, zSlots};
            for (std::size_t y = 0; y < ySlots; ++y)
            {
                for (std::size_t z = 0; z < zSlots; ++z)
                {
                    const auto index = transverseIndex[y * zSlots + z];
                    if (index != missing)
                    {
                        waves[0][y] = transverseWaves[index][0];
                        waves[1][z] = transverseWaves[index][1];
                    }
                }
            }
            for (std::size_t axis = 0; axis < 2; ++axis)
            {
                auto &values = coordinates[axis];
                values.reserve(m_points.size());
                for (const auto &point : m_points)
                {
                    values.push_back(point[axis + 1]);
                }
                std::sort(values.begin(), values.end());
                values.erase(std::unique(values.begin(), values.end()),
                             values.end());
                uniqueCoordinates[axis] = values.size();
                std::size_t count = 0;
                if (!Product(values.size(), slots[axis], count) ||
                    !AddPayload(count, sizeof(Complex), maxCacheBytes, bytes))
                {
                    m_method = "ReferenceCacheLimit";
                    return;
                }
                auto &phases = axisPhases[axis];
                phases.resize(count);
                if (phases.capacity() > phases.size() &&
                    !AddPayload(phases.capacity() - phases.size(),
                                sizeof(Complex), maxCacheBytes, bytes))
                {
                    m_method = "ReferenceCacheLimit";
                    return;
                }
                for (std::size_t coordinate = 0; coordinate < values.size();
                     ++coordinate)
                {
                    const long double displacement =
                        static_cast<long double>(values[coordinate]) -
                        config.origin[axis + 1];
                    for (std::size_t wave = 0; wave < slots[axis]; ++wave)
                    {
                        phases[coordinate * slots[axis] + wave] =
                            Phase(waves[axis][wave] * displacement);
                    }
                }
                for (std::size_t point = 0; point < m_points.size(); ++point)
                {
                    pointIndex[point][axis] = static_cast<std::size_t>(
                        std::lower_bound(values.begin(), values.end(),
                                         m_points[point][axis + 1]) -
                        values.begin());
                }
            }
            std::array<std::size_t, 2> workspace;
            if (!Product(uniqueCoordinates[0], zSlots, workspace[0]) ||
                !Product(uniqueCoordinates[1], ySlots, workspace[1]))
            {
                m_method = "ReferenceCacheLimit";
                return;
            }
            firstAxis = workspace[0] <= workspace[1] ? 0 : 1;
            if (!AddPayload(workspace[firstAxis], sizeof(HITComplexVector),
                            maxCacheBytes, bytes))
            {
                m_method = "ReferenceCacheLimit";
                return;
            }
        }
        m_xIndex = std::move(xIndex);
        m_transverseIndex = std::move(transverseIndex);
        m_transverseOffset = {maxIndex[1], maxIndex[2]};
        m_transverseStride = zSlots;
        m_xWaves = std::move(xWaves);
        m_transverseWaves = std::move(transverseWaves);
        m_spatialPhases = std::move(spatialPhases);
        m_axisPhases = std::move(axisPhases);
        m_pointIndex = std::move(pointIndex);
        m_uniqueCoordinates = uniqueCoordinates;
        m_firstAxis = firstAxis;
        m_cacheBytes = bytes;
        m_separable = !fullCache;
        m_grouped = true;
        m_method = fullCache ? "GroupedTransverse" : "GroupedSeparable";
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

std::size_t HITPlaneEvaluator::TransverseBin(int y, int z) const
{
    const auto iy = static_cast<std::size_t>(
        static_cast<std::int64_t>(y) +
        static_cast<std::int64_t>(m_transverseOffset[0]));
    const auto iz = static_cast<std::size_t>(
        static_cast<std::int64_t>(z) +
        static_cast<std::int64_t>(m_transverseOffset[1]));
    return m_transverseIndex[iy * m_transverseStride + iz];
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

void HITPlaneEvaluator::EvaluateSeparable(
    const std::vector<HITComplexVector> &amplitude,
    std::vector<HITVector> &velocities) const
{
    const std::array<std::size_t, 2> slots{
        2 * m_transverseOffset[0] + 1, m_transverseStride};
    const auto first = m_firstAxis;
    const auto second = 1 - first;
    const auto firstSlots = slots[first];
    const auto secondSlots = slots[second];
    // Constructor preflight bounds this workspace and picks the direction
    // requiring fewer rows. No Cartesian-product assumption is made about
    // the supplied planar points.
    std::vector<HITComplexVector> intermediate(
        m_uniqueCoordinates[first] * secondSlots);
    const auto missing = std::numeric_limits<std::size_t>::max();
    for (std::size_t coordinate = 0;
         coordinate < m_uniqueCoordinates[first]; ++coordinate)
    {
        for (std::size_t other = 0; other < secondSlots; ++other)
        {
            auto &sum = intermediate[coordinate * secondSlots + other];
            for (std::size_t wave = 0; wave < firstSlots; ++wave)
            {
                const auto index = first == 0 ?
                    wave * m_transverseStride + other :
                    other * m_transverseStride + wave;
                const auto bin = m_transverseIndex[index];
                if (bin == missing)
                {
                    continue;
                }
                const auto phase =
                    m_axisPhases[first][coordinate * firstSlots + wave];
                for (int component = 0; component < 3; ++component)
                {
                    sum[component] += amplitude[bin][component] * phase;
                }
            }
        }
    }
    const auto speed = m_generator.GetConfig().meanSpeed;
    for (std::size_t point = 0; point < m_points.size(); ++point)
    {
        std::array<long double, 3> sum{speed, 0.0L, 0.0L};
        const auto row = m_pointIndex[point][first] * secondSlots;
        const auto phases = m_pointIndex[point][second] * secondSlots;
        for (std::size_t wave = 0; wave < secondSlots; ++wave)
        {
            const auto phase = m_axisPhases[second][phases + wave];
            for (int component = 0; component < 3; ++component)
            {
                sum[component] +=
                    (intermediate[row + wave][component] * phase).real();
            }
        }
        for (int component = 0; component < 3; ++component)
        {
            velocities[point][component] = Velocity(sum[component]);
        }
    }
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
    const auto &modes = m_generator.GetModes();
    velocities.assign(m_points.size(), HITVector{config.meanSpeed, 0.0, 0.0});
    if (m_points.empty() || modes.empty())
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
    for (const auto &mode : modes)
    {
        const auto xPhase = m_xIndex[static_cast<std::size_t>(mode.index[0])];
        const auto positiveBin = TransverseBin(mode.index[1], mode.index[2]);
        const auto negativeBin = TransverseBin(-mode.index[1], -mode.index[2]);
        for (int component = 0; component < 3; ++component)
        {
            const auto value = mode.coefficient[component] * xPhases[xPhase];
            if (positiveBin == negativeBin)
            {
                amplitude[positiveBin][component] +=
                    Complex(2.0 * value.real(), 0.0);
            }
            else
            {
                amplitude[positiveBin][component] += value;
                amplitude[negativeBin][component] += std::conj(value);
            }
        }
    }
    if (m_separable)
    {
        EvaluateSeparable(amplitude, velocities);
        return;
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
