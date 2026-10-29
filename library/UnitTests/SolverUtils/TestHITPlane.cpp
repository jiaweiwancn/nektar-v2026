///////////////////////////////////////////////////////////////////////////////
// File: TestHITPlane.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Turbulence/HITPlaneEvaluator.h>
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace Nektar::SolverUtils::HITPlaneTests
{
namespace
{
const NekDouble twoPi = 2.0 * std::acos(-1.0);

HITGeneratorConfig Configuration(NekDouble length = 1.0,
                                 NekDouble velocity = 1.0)
{
    HITGeneratorConfig config;
    config.lengths = {twoPi * length, twoPi * length / 1.3,
                      twoPi * length / 0.9};
    config.origin = {0.7 * length, -0.4 * length, 1.2 * length};
    config.counts = {10, 10, 10};
    config.seed = 20261028;
    config.meanSpeed = 1.3 * velocity;
    config.kMax = 3.3 / length;
    config.taperStart = 2.7 / length;
    config.normalization = HITNormalization::TargetRMS;
    config.targetRMS = 0.05 * velocity;
    config.timeStep = 0.001 * length / velocity;
    config.recordDuration = 0.5 * length / velocity;
    return config;
}

std::vector<HITVector> Points(NekDouble length = 1.0)
{
    std::vector<HITVector> result;
    for (int y = 0; y < 5; ++y)
    {
        for (int z = 0; z < 8; ++z)
        {
            result.push_back({0.21 * length, (-0.43 + 0.37 * y) * length,
                              (0.67 + z * twoPi / 8.0) * length});
        }
    }
    result.push_back(result.front()); // Duplicate element-interface point.
    return result;
}

std::vector<HITVector> Compare(const HITGenerator &generator,
                              const HITPlaneEvaluator &plane,
                              const std::vector<HITVector> &points,
                              NekDouble time, bool exact = false)
{
    std::vector<HITVector> actual;
    plane.Evaluate(time, actual);
    BOOST_REQUIRE_EQUAL(actual.size(), points.size());
    const auto &config = generator.GetConfig();
    for (std::size_t point = 0; point < points.size(); ++point)
    {
        const auto expected = generator.Evaluate(points[point], time).velocity;
        for (int component = 0; component < 3; ++component)
        {
            if (exact)
            {
                BOOST_CHECK_EQUAL(actual[point][component], expected[component]);
            }
            else
            {
                const auto scale = std::max(
                    component == 0 ? config.meanSpeed : config.targetRMS,
                    std::abs(expected[component]));
                BOOST_REQUIRE(scale > 0.0);
                BOOST_CHECK_LE(std::abs(actual[point][component] -
                                        expected[component]) / scale,
                               5e-12);
            }
        }
    }
    return actual;
}
} // namespace

BOOST_AUTO_TEST_SUITE(HITCore)

BOOST_AUTO_TEST_CASE(PlaneTransverseConjugateGroupsAndRepeatedTimes)
{
    const auto config = Configuration();
    const HITGenerator generator(config, HITSpectrum::VonKarman(0.8, 1.0));
    bool oblique = false, negativeY = false, negativeZ = false;
    bool zeroX = false, zeroTransverse = false;
    for (const auto &mode : generator.GetModes())
    {
        oblique |= mode.index[0] && mode.index[1] && mode.index[2];
        negativeY |= mode.index[1] < 0;
        negativeZ |= mode.index[2] < 0;
        zeroX |= mode.index[0] == 0;
        zeroTransverse |= mode.index[1] == 0 && mode.index[2] == 0;
    }
    BOOST_REQUIRE(oblique && negativeY && negativeZ && zeroX && zeroTransverse);
    const auto points = Points();
    const HITPlaneEvaluator plane(generator, points);
    BOOST_REQUIRE(plane.IsGrouped());
    BOOST_CHECK_EQUAL(plane.Method(), "GroupedTransverse");
    BOOST_CHECK_GT(plane.CacheBytes(), 0U);
    BOOST_CHECK_LE(plane.CacheBytes(), 64U * 1024U * 1024U);
    const auto first = Compare(generator, plane, points, 0.0);
    const auto middle = Compare(generator, plane, points, 0.2);
    Compare(generator, plane, points, config.recordDuration);
    BOOST_CHECK(Compare(generator, plane, points, 0.2) == middle);
    BOOST_CHECK(Compare(generator, plane, points, 0.0) == first);
    BOOST_CHECK(first.front() == first.back());
}

BOOST_AUTO_TEST_CASE(PlaneZeroRMSEmptyAndFiniteCoordinates)
{
    auto config = Configuration();
    config.targetRMS = 0.0;
    const HITGenerator generator(config, HITSpectrum::VonKarman(0.8, 1.0));
    const auto points = Points();
    const HITPlaneEvaluator plane(generator, points, 0);
    BOOST_REQUIRE(plane.IsGrouped());
    BOOST_CHECK_EQUAL(plane.Method(), "MeanOnly");
    BOOST_CHECK_EQUAL(plane.CacheBytes(), 0U);
    Compare(generator, plane, points, 0.0, true);
    Compare(generator, plane, points, config.recordDuration, true);
    const HITPlaneEvaluator empty(generator, {});
    std::vector<HITVector> values(3);
    empty.Evaluate(0.0, values);
    BOOST_CHECK(values.empty());
    BOOST_CHECK_THROW(empty.Evaluate(1.0, values), std::invalid_argument);
    auto invalid = points;
    invalid[1][1] = std::numeric_limits<NekDouble>::quiet_NaN();
    BOOST_CHECK_THROW(HITPlaneEvaluator(generator, invalid),
                      std::invalid_argument);
    invalid[1][1] = 0.0;
    invalid[1][0] = std::numeric_limits<NekDouble>::infinity();
    BOOST_CHECK_THROW(HITPlaneEvaluator(generator, invalid),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(PlaneNonplanarAndCacheLimitUseReference)
{
    const auto config = Configuration();
    const HITGenerator generator(config, HITSpectrum::VonKarman(0.8, 1.0));
    auto points = Points();
    const HITPlaneEvaluator grouped(generator, points);
    BOOST_REQUIRE(grouped.IsGrouped());
    const HITPlaneEvaluator tooSmall(generator, points, grouped.CacheBytes() - 1);
    BOOST_CHECK(!tooSmall.IsGrouped());
    BOOST_CHECK_EQUAL(tooSmall.CacheBytes(), 0U);
    BOOST_CHECK_EQUAL(tooSmall.Method(), "ReferenceCacheLimit");
    Compare(generator, tooSmall, points, 0.2, true);
    const HITPlaneEvaluator exactBudget(generator, points, grouped.CacheBytes());
    BOOST_CHECK(exactBudget.IsGrouped());
    Compare(generator, exactBudget, points, 0.2);
    const HITPlaneEvaluator noCache(generator, points, 0);
    BOOST_CHECK(!noCache.IsGrouped());
    Compare(generator, noCache, points, 0.2, true);

    points.back()[0] += 0.01;
    const HITPlaneEvaluator nonplanar(generator, points);
    BOOST_CHECK(!nonplanar.IsGrouped());
    BOOST_CHECK_EQUAL(nonplanar.Method(), "ReferenceNonplanar");
    BOOST_CHECK_EQUAL(nonplanar.CacheBytes(), 0U);
    Compare(generator, nonplanar, points, 0.0, true);
    Compare(generator, nonplanar, points, 0.2, true);
}

BOOST_AUTO_TEST_CASE(PlaneSharesRecordValidationWithLargeTimeOrigin)
{
    auto config = Configuration();
    config.timeOrigin = 1e12;
    const HITGenerator generator(config, HITSpectrum::VonKarman(0.8, 1.0));
    const auto points = Points();
    const HITPlaneEvaluator plane(generator, points);
    BOOST_REQUIRE(plane.IsGrouped());
    for (auto offset : {0.0, 0.125, 0.5, 0.125})
    {
        Compare(generator, plane, points, config.timeOrigin + offset);
    }
    std::vector<HITVector> values;
    for (auto time : {config.timeOrigin - 0.01,
                      config.timeOrigin + config.recordDuration + 0.01,
                      std::numeric_limits<NekDouble>::quiet_NaN(),
                      std::numeric_limits<NekDouble>::infinity()})
    {
        BOOST_CHECK_THROW(plane.Evaluate(time, values), std::invalid_argument);
    }
}

BOOST_AUTO_TEST_CASE(PlaneLargeSpatialPhasesUseReference)
{
    const auto config = Configuration();
    const HITGenerator generator(config, HITSpectrum::VonKarman(0.8, 1.0));
    for (int direction : {0, 1, 2})
    {
        auto points = Points();
        for (auto &point : points)
        {
            point[direction] += 1e15;
        }
        const HITPlaneEvaluator plane(generator, points);
        BOOST_CHECK(!plane.IsGrouped());
        BOOST_CHECK_EQUAL(plane.Method(), "ReferenceLargePhase");
        Compare(generator, plane, points, 0.0, true);
        Compare(generator, plane, points, 0.2, true);
    }
    // Safety is determined over the declared record, not only its first call.
    auto longRecord = config;
    longRecord.recordDuration = 1e12;
    longRecord.allowRepeat = true;
    const HITGenerator longGenerator(longRecord,
                                      HITSpectrum::VonKarman(0.8, 1.0));
    const auto points = Points();
    const HITPlaneEvaluator longPlane(longGenerator, points);
    BOOST_CHECK(!longPlane.IsGrouped());
    BOOST_CHECK_EQUAL(longPlane.Method(), "ReferenceLargePhase");
    Compare(longGenerator, longPlane, points, 0.0, true);
    Compare(longGenerator, longPlane, points, 1e11, true);
}

BOOST_AUTO_TEST_CASE(PlaneScaledUnitsAndTabulatedSpectrum)
{
    for (NekDouble length : {1e-6, 1.0, 1e6})
    {
        for (NekDouble velocity : {1e-8, 1.0, 1e8})
        {
            BOOST_TEST_CONTEXT("length=" << length << ", velocity=" << velocity)
            {
                const auto config = Configuration(length, velocity);
                const HITGenerator generator(config,
                    HITSpectrum::VonKarman(0.8 * length, 1.0));
                const auto points = Points(length);
                const HITPlaneEvaluator plane(generator, points);
                BOOST_REQUIRE(plane.IsGrouped());
                Compare(generator, plane, points, 0.23 * length / velocity);
                Compare(generator, plane, points, config.recordDuration);
            }
        }
    }
    std::istringstream input(
        "# NEKTAR_HIT_SPECTRUM_V1 KIND=Shell3D K_UNITS=rad/m E_UNITS=m3/s2\n"
        "0.5 0.02\n1 0.1\n1.5 0.2\n2 0\n");
    const auto spectrum = HITSpectrum::ReadTable(
        input, HITInterpolation::Linear, HITLowWavenumber::Zero,
        HITHighWavenumber::Zero);
    const HITGenerator generator(Configuration(), spectrum);
    const auto points = Points();
    const HITPlaneEvaluator plane(generator, points);
    BOOST_REQUIRE(plane.IsGrouped());
    Compare(generator, plane, points, 0.2);
    Compare(generator, plane, points, 0.0);
}

BOOST_AUTO_TEST_CASE(PlaneDerivativeRangePreservesReferenceExceptions)
{
    struct Scenario
    {
        const char *name;
        NekDouble length, speed, rms, dt, duration;
        bool partialProductOnly;
    };
    const std::array<Scenario, 2> scenarios{{
        {"unrepresentable derivative", 1.0, 1e250, 1e100, 1e-252,
         1e-251, false},
        {"overflowing intermediate product", 1e10, 1e308, 1e-100, 1e-307,
         1e-306, true}}};
    for (const auto &scenario : scenarios)
    {
        BOOST_TEST_CONTEXT(scenario.name)
        {
            auto config = Configuration(scenario.length);
            config.meanSpeed = scenario.speed;
            config.targetRMS = scenario.rms;
            config.timeStep = scenario.dt;
            config.recordDuration = scenario.duration;
            const HITGenerator generator(config,
                HITSpectrum::VonKarman(0.8 * scenario.length, 1.0));
            const std::vector<HITVector> points{config.origin};
            // Every phase is zero here. Velocities remain finite even when
            // the reference's left-associated time-derivative product fails.
            HITVector velocity{config.meanSpeed, 0.0, 0.0};
            long double maxWave = 0.0L, maxCoefficient = 0.0L;
            for (const auto &mode : generator.GetModes())
            {
                for (int component = 0; component < 3; ++component)
                {
                    velocity[component] +=
                        2.0 * mode.coefficient[component].real();
                    maxWave = std::max(maxWave,
                        std::abs(static_cast<long double>(
                            mode.wavevector[component])));
                    maxCoefficient = std::max(maxCoefficient,
                        static_cast<long double>(
                            std::abs(mode.coefficient[component])));
                }
            }
            for (const auto value : velocity)
            {
                BOOST_REQUIRE(std::isfinite(value));
            }
            if (scenario.partialProductOnly)
            {
                // Isolate the partial-product guard: the complete paired
                // derivative bound fits, but 2*U overflows before the small
                // wavenumber and coefficient are multiplied. Timesteps and
                // record lengths stay in the normal floating-point range.
                const long double bound =
                    (2.0L * generator.GetModes().size() * maxWave *
                     maxCoefficient) * config.meanSpeed;
                BOOST_REQUIRE(std::isfinite(bound));
                BOOST_CHECK_LE(bound,
                    static_cast<long double>(
                        std::numeric_limits<NekDouble>::max()) / 16.0L);
                BOOST_CHECK_GT(2.0L * config.meanSpeed,
                    static_cast<long double>(
                        std::numeric_limits<NekDouble>::max()));
            }
            BOOST_CHECK_THROW(generator.Evaluate(config.origin, 0.0),
                              std::overflow_error);
            const HITPlaneEvaluator unlimited(
                generator, points, std::numeric_limits<std::size_t>::max());
            const HITPlaneEvaluator forcedFallback(generator, points, 0);
            std::vector<HITVector> values;
            for (const auto *plane : {&unlimited, &forcedFallback})
            {
                BOOST_CHECK(!plane->IsGrouped());
                BOOST_CHECK_EQUAL(plane->Method(), "ReferenceDerivativeRange");
                BOOST_CHECK_EQUAL(plane->CacheBytes(), 0U);
                BOOST_CHECK_THROW(plane->Evaluate(0.0, values),
                                  std::overflow_error);
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
} // namespace Nektar::SolverUtils::HITPlaneTests
