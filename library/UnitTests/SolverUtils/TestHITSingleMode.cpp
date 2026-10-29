///////////////////////////////////////////////////////////////////////////////
// File: TestHITSingleMode.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Turbulence/HITGenerator.h>
#include <boost/test/unit_test.hpp>

#include <array>
#include <cmath>
#include <complex>
#include <limits>

namespace Nektar::SolverUtils::HITInvariantTests
{
namespace
{
const NekDouble twoPi = 2.0 * std::acos(-1.0);

long double Norm(const HITVector &values)
{
    long double sum = 0.0L;
    for (auto value : values)
    {
        sum += static_cast<long double>(value) * value;
    }
    return std::sqrt(sum);
}

long double GradientNorm(const HITTensor &gradient)
{
    long double sum = 0.0L;
    for (const auto &row : gradient)
    {
        for (auto value : row)
        {
            sum += static_cast<long double>(value) * value;
        }
    }
    return std::sqrt(sum);
}

long double NormalizedResidual(long double numerator, long double denominator)
{
    if (denominator == 0.0L)
    {
        return numerator == 0.0L
                   ? 0.0L
                   : std::numeric_limits<long double>::infinity();
    }
    return std::abs(numerator) / denominator;
}

HITGeneratorConfig ScaledConfiguration(NekDouble length, NekDouble velocity)
{
    HITGeneratorConfig config;
    config.lengths = {twoPi * length, twoPi * length / 1.5,
                      twoPi * length / 1.2};
    config.origin = {0.7 * length, -0.4 * length, 1.2 * length};
    config.counts = {8, 8, 8};
    config.seed = 20261028;
    config.meanSpeed = 1.3 * velocity;
    config.kMax = 3.3 / length;
    config.taperStart = 2.8 / length;
    config.normalization = HITNormalization::TargetRMS;
    config.targetRMS = 0.05 * velocity;
    config.timeStep = 0.001 * length / velocity;
    config.recordDuration = 0.5 * length / velocity;
    return config;
}
} // namespace

BOOST_AUTO_TEST_SUITE(HITCore)

BOOST_AUTO_TEST_CASE(DimensionlessSolenoidalityAcrossPhysicalScales)
{
    // No dimensionful max(1,scale) floor: reducing either velocity or physical
    // length must not make a longitudinal component easier to hide.
    for (NekDouble length : {1e-6, 1.0, 1e6})
    {
        for (NekDouble velocity : {1e-8, 1.0, 1e8})
        {
            BOOST_TEST_CONTEXT("length=" << length << ", velocity=" << velocity)
            {
                const auto config = ScaledConfiguration(length, velocity);
                const HITGenerator generator(config,
                    HITSpectrum::VonKarman(0.8 * length, 1.0));
                BOOST_REQUIRE(!generator.GetModes().empty());
                for (const auto &mode : generator.GetModes())
                {
                    std::complex<long double> dot{};
                    long double coefficientSquared = 0.0L;
                    for (int direction = 0; direction < 3; ++direction)
                    {
                        const auto coefficient = mode.coefficient[direction];
                        const std::complex<long double> wide(
                            coefficient.real(), coefficient.imag());
                        dot += static_cast<long double>(mode.wavevector[direction]) * wide;
                        coefficientSquared += std::norm(wide);
                    }
                    const auto scale = Norm(mode.wavevector) *
                                       std::sqrt(coefficientSquared);
                    BOOST_CHECK_LE(NormalizedResidual(std::abs(dot), scale),
                                   1e-12L);
                }
                const auto sample = generator.Evaluate(
                    {0.21 * length, -0.43 * length, 0.67 * length},
                    0.2 * length / velocity);
                const long double trace =
                    static_cast<long double>(sample.gradient[0][0]) +
                    sample.gradient[1][1] + sample.gradient[2][2];
                const auto scale = GradientNorm(sample.gradient);
                BOOST_REQUIRE(std::isfinite(scale) && scale > 0.0L);
                BOOST_CHECK_LE(NormalizedResidual(trace, scale), 1e-12L);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(DimensionlessAnalyticDerivativesAcrossPhysicalScales)
{
    for (NekDouble length : {1e-6, 1.0, 1e6})
    {
        for (NekDouble velocity : {1e-8, 1.0, 1e8})
        {
            BOOST_TEST_CONTEXT("length=" << length << ", velocity=" << velocity)
            {
                const auto config = ScaledConfiguration(length, velocity);
                const HITGenerator generator(config,
                    HITSpectrum::VonKarman(0.8 * length, 1.0));
                const HITVector point{
                    0.21 * length, -0.43 * length, 0.67 * length};
                const auto timeScale = length / velocity;
                const NekDouble time = 0.2 * timeScale;
                const auto sample = generator.Evaluate(point, time);
                const NekDouble h = 1e-5 * length;
                long double spatialErrorSquared = 0.0L;
                for (int direction = 0; direction < 3; ++direction)
                {
                    auto plus = point, minus = point;
                    plus[direction] += h;
                    minus[direction] -= h;
                    const auto forward = generator.Evaluate(plus, time);
                    const auto backward = generator.Evaluate(minus, time);
                    for (int component = 0; component < 3; ++component)
                    {
                        const long double difference =
                            (forward.velocity[component] - backward.velocity[component]) /
                                (2.0 * h) - sample.gradient[component][direction];
                        spatialErrorSquared += difference * difference;
                    }
                }
                const auto gradientScale = GradientNorm(sample.gradient);
                BOOST_REQUIRE(gradientScale > 0.0L);
                BOOST_CHECK_LE(NormalizedResidual(std::sqrt(spatialErrorSquared),
                                                  gradientScale), 1e-8L);

                const NekDouble timeStep = 1e-5 * timeScale;
                const auto forward = generator.Evaluate(point, time + timeStep);
                const auto backward = generator.Evaluate(point, time - timeStep);
                HITVector timeError{}, convectiveResidual{}, xGradient{};
                for (int component = 0; component < 3; ++component)
                {
                    timeError[component] =
                        (forward.velocity[component] - backward.velocity[component]) /
                            (2.0 * timeStep) - sample.timeDerivative[component];
                    xGradient[component] = sample.gradient[component][0];
                    convectiveResidual[component] = sample.timeDerivative[component] +
                        config.meanSpeed * xGradient[component];
                }
                const auto temporalScale = Norm(sample.timeDerivative);
                BOOST_REQUIRE(temporalScale > 0.0L);
                BOOST_CHECK_LE(NormalizedResidual(Norm(timeError), temporalScale),
                               1e-8L);
                BOOST_CHECK_LE(NormalizedResidual(
                    Norm(convectiveResidual), temporalScale +
                        config.meanSpeed * Norm(xGradient)), 1e-12L);

                auto followedPoint = point;
                const NekDouble shift = 0.07 * timeScale;
                followedPoint[0] += config.meanSpeed * shift;
                const auto followed = generator.Evaluate(followedPoint, time + shift);
                HITVector difference{};
                for (int component = 0; component < 3; ++component)
                {
                    difference[component] = followed.velocity[component] -
                                            sample.velocity[component];
                }
                BOOST_CHECK_LE(NormalizedResidual(Norm(difference), config.targetRMS),
                               1e-12L);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(ZeroFieldHasExplicitZeroInvariantNormalization)
{
    auto config = ScaledConfiguration(1e6, 1e-8);
    config.targetRMS = 0.0;
    const HITGenerator generator(config,
        HITSpectrum::VonKarman(0.8e6, 1.0));
    BOOST_CHECK(generator.GetModes().empty());
    const auto sample = generator.Evaluate({0.21e6, -0.43e6, 0.67e6},
                                           0.2e6 / 1e-8);
    const auto scale = GradientNorm(sample.gradient);
    const long double trace = static_cast<long double>(sample.gradient[0][0]) +
                              sample.gradient[1][1] + sample.gradient[2][2];
    BOOST_CHECK_EQUAL(scale, 0.0L);
    BOOST_CHECK_EQUAL(trace, 0.0L);
    BOOST_CHECK_EQUAL(NormalizedResidual(trace, scale), 0.0L);
    BOOST_CHECK_EQUAL(Norm(sample.timeDerivative), 0.0L);
    BOOST_CHECK_EQUAL(sample.velocity[0], config.meanSpeed);
    BOOST_CHECK_EQUAL(sample.velocity[1], 0.0);
    BOOST_CHECK_EQUAL(sample.velocity[2], 0.0);
}

BOOST_AUTO_TEST_CASE(SeededSingleTransversePairMatchesKnownWave)
{
    HITGeneratorConfig config;
    config.lengths = {twoPi, 1.0, 1.0};
    config.counts = {4, 4, 4};
    config.seed = 20261028;
    config.kMax = 1.55;
    config.taperStart = 1.5;
    config.normalization = HITNormalization::TargetRMS;
    config.targetRMS = 0.05;
    config.timeStep = 0.005;
    config.recordDuration = 0.1;
    const HITGenerator generator(config, HITSpectrum::VonKarman(0.8, 1.0));
    BOOST_REQUIRE_EQUAL(generator.GetModes().size(), 1U);
    const auto &mode = generator.GetModes()[0];
    BOOST_CHECK((mode.index == HITIndex{1, 0, 0}));
    BOOST_CHECK_EQUAL(HITGenerator::ModeId(mode.index, config.counts), 58U);
    BOOST_CHECK_EQUAL(mode.wavevector[0], 1.0);
    BOOST_CHECK_EQUAL(mode.wavevector[1], 0.0);
    BOOST_CHECK_EQUAL(mode.wavevector[2], 0.0);
    // Independent SplitMix64/Box-Muller oracle in generate_m5_single_mode.py.
    const HITComplexVector expected{{{0.0, 0.0},
        {-0.011280905803074937, -0.0015420401494205553},
        {-0.017914370023552, -0.030473381145025093}}};
    for (int component = 0; component < 3; ++component)
    {
        BOOST_CHECK_LE(std::abs(mode.coefficient[component] - expected[component]) /
                           config.targetRMS, 2e-13);
    }
    const auto sample = generator.Evaluate({0.37, 0.37, 0.2}, 0.073);
    const HITVector expectedVelocity{
        1.0, -0.0206714623224743, -0.016423867998451948};
    const HITVector expectedXGradient{
        0.0, 0.009551834362656447, 0.06876382047975711};
    for (int component = 0; component < 3; ++component)
    {
        BOOST_CHECK_LE(std::abs(sample.velocity[component] -
                                expectedVelocity[component]) / config.targetRMS,
                       3e-13);
        BOOST_CHECK_LE(std::abs(sample.gradient[component][0] -
                                expectedXGradient[component]) / config.targetRMS,
                       3e-13);
        BOOST_CHECK_LE(std::abs(sample.timeDerivative[component] +
                                expectedXGradient[component]) / config.targetRMS,
                       3e-13);
        BOOST_CHECK_EQUAL(sample.gradient[component][1], 0.0);
        BOOST_CHECK_EQUAL(sample.gradient[component][2], 0.0);
    }
}

BOOST_AUTO_TEST_SUITE_END()
} // namespace Nektar::SolverUtils::HITInvariantTests
