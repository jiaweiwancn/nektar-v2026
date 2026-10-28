///////////////////////////////////////////////////////////////////////////////
// File: TestHIT.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Turbulence/HITGenerator.h>
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <fstream>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace Nektar::SolverUtils::HITUnitTests
{
namespace
{
constexpr NekDouble pi = 3.141592653589793238462643383279502884;
const std::string header =
    "# NEKTAR_HIT_SPECTRUM_V1 KIND=Shell3D K_UNITS=rad/m E_UNITS=m3/s2\n";

void Near(NekDouble actual, NekDouble expected, NekDouble tolerance = 2e-13)
{
    BOOST_CHECK_SMALL(actual - expected,
                      tolerance * std::max(NekDouble(1.0), std::abs(expected)));
}

NekDouble Taper(NekDouble k, const HITGeneratorConfig &config)
{
    if (k >= config.kMax)
    {
        return 0.0;
    }
    if (k <= config.taperStart)
    {
        return 1.0;
    }
    return 0.5 * (1.0 + std::cos(pi * (k - config.taperStart) /
                               (config.kMax - config.taperStart)));
}

struct ReferenceMoments
{
    HITTensor covariance{};
    NekDouble energy = 0.0;
    NekDouble dissipation = 0.0;
    NekDouble energyVariance = 0.0;
    NekDouble fourthCumulant = 0.0;
    std::size_t fullModes = 0;
};

// Integrate the Batchelor tensor over the FULL lattice, independently of the
// implementation's half-lattice generation and cross-product factorization.
ReferenceMoments TensorMoments(const HITGeneratorConfig &config,
                              const HITSpectrum &spectrum,
                              NekDouble viscosity, NekDouble energyScale = 1.0)
{
    ReferenceMoments result;
    const NekDouble deltaVolume = std::pow(2.0 * pi, 3) /
                                  (config.lengths[0] * config.lengths[1] *
                                   config.lengths[2]);
    for (int x = -config.counts[0] / 2 + 1; x < config.counts[0] / 2; ++x)
    {
        for (int y = -config.counts[1] / 2 + 1; y < config.counts[1] / 2; ++y)
        {
            for (int z = -config.counts[2] / 2 + 1;
                 z < config.counts[2] / 2; ++z)
            {
                HITVector k{2.0 * pi * x / config.lengths[0],
                            2.0 * pi * y / config.lengths[1],
                            2.0 * pi * z / config.lengths[2]};
                NekDouble k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
                if (k2 == 0.0 || std::sqrt(k2) >= config.kMax)
                {
                    continue;
                }
                NekDouble lambda = energyScale * deltaVolume *
                                   spectrum.Evaluate(std::sqrt(k2)) *
                                   Taper(std::sqrt(k2), config) /
                                   (4.0 * pi * k2);
                if (lambda == 0.0)
                {
                    continue;
                }
                ++result.fullModes;
                for (int i = 0; i < 3; ++i)
                {
                    for (int j = 0; j < 3; ++j)
                    {
                        result.covariance[i][j] +=
                            lambda * ((i == j ? 1.0 : 0.0) - k[i] * k[j] / k2);
                    }
                }
                result.energy += lambda;
                result.dissipation += 2.0 * viscosity * k2 * lambda;
                // Each independent proper complex transverse coefficient has
                // two equal eigenvalues lambda. Conjugate partners count once
                // in K. Thus Var(K)=sum_full lambda^2 and kappa4(K)=sum_full
                // 6*lambda^4, without estimating either from the realization.
                result.energyVariance += lambda * lambda;
                result.fourthCumulant += 6.0 * std::pow(lambda, 4);
            }
        }
    }
    return result;
}

void CheckSameSample(const HITSample &left, const HITSample &right)
{
    for (int i = 0; i < 3; ++i)
    {
        BOOST_CHECK_EQUAL(left.velocity[i], right.velocity[i]);
        BOOST_CHECK_EQUAL(left.timeDerivative[i], right.timeDerivative[i]);
        for (int j = 0; j < 3; ++j)
        {
            BOOST_CHECK_EQUAL(left.gradient[i][j], right.gradient[i][j]);
        }
    }
}
} // namespace

BOOST_AUTO_TEST_SUITE(HITCore)

BOOST_AUTO_TEST_CASE(RNGKnownAnswersAndOpenUniformInterval)
{
    BOOST_CHECK_EQUAL(std::string(HITGenerator::Algorithm()),
                      "SplitMix64BoxMuller52V1");
    BOOST_CHECK_EQUAL(HITGenerator::Mix64(0), UINT64_C(0xe220a8397b1dcdaf));
    BOOST_CHECK_EQUAL(HITGenerator::Mix64(1), UINT64_C(0x910a2dec89025cc1));
    BOOST_CHECK_EQUAL(HITGenerator::Mix64(2), UINT64_C(0x975835de1c9756ce));
    BOOST_CHECK_EQUAL(HITGenerator::Uniform52(0), std::ldexp(1.0, -53));
    BOOST_CHECK_EQUAL(HITGenerator::Uniform52(UINT64_MAX),
                      1.0 - std::ldexp(1.0, -53));
    BOOST_CHECK(HITGenerator::Uniform52(0) > 0.0);
    BOOST_CHECK(HITGenerator::Uniform52(UINT64_MAX) < 1.0);
    BOOST_CHECK_EQUAL(HITGenerator::ModeId({0, 0, 1}, {8, 8, 8}), 293U);
    BOOST_CHECK_EQUAL(HITGenerator::ModeId({1, -2, 3}, {8, 8, 8}), 343U);
    const auto words = HITGenerator::CounterWords(UINT64_MAX, 0);
    // Adding q wraps the seed modulo 2^64, rather than saturating it.
    BOOST_CHECK_EQUAL(words[1], UINT64_C(0xe220a8397b1dcdaf));
    BOOST_CHECK_EQUAL(words[2], UINT64_C(0x910a2dec89025cc1));
    BOOST_CHECK_EQUAL(words[3], UINT64_C(0x975835de1c9756ce));
    BOOST_CHECK_THROW(HITGenerator::CounterWords(0, UINT64_MAX),
                      std::overflow_error);
}

BOOST_AUTO_TEST_CASE(StrictSeedAndModeIndexValidation)
{
    BOOST_CHECK_EQUAL(HITGenerator::ParseSeed("0"), 0U);
    BOOST_CHECK_EQUAL(HITGenerator::ParseSeed("18446744073709551615"), UINT64_MAX);
    for (const std::string text : {"", "-1", "+1", "1.0", " 1", "1 ", "1e3", "one"})
    {
        BOOST_CHECK_THROW(HITGenerator::ParseSeed(text), std::invalid_argument);
    }
    BOOST_CHECK_THROW(HITGenerator::ParseSeed("18446744073709551616"),
                      std::overflow_error);
    BOOST_CHECK_THROW(HITGenerator::ModeId({4, 0, 0}, {8, 8, 8}),
                      std::invalid_argument);
    BOOST_CHECK_THROW(HITGenerator::ModeId({1, 0, 0}, {7, 8, 8}),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(VonKarmanSpectrumFormulaAndRange)
{
    const NekDouble length = 0.3, amplitude = 1.7;
    auto spectrum = HITSpectrum::VonKarman(length, amplitude);
    for (NekDouble k : {0.1, 1.0, 5.0, 100.0})
    {
        const NekDouble lk = length * k;
        const NekDouble expected = amplitude * std::pow(length, 5.0 / 3.0) *
                                   std::pow(lk, 4) /
                                   std::pow(1.0 + lk * lk, 17.0 / 6.0);
        BOOST_CHECK_CLOSE_FRACTION(spectrum.Evaluate(k), expected, 5e-14);
    }
    BOOST_CHECK_EQUAL(spectrum.Evaluate(0.0), 0.0);
    BOOST_CHECK(std::isfinite(spectrum.Evaluate(1e100)));
    BOOST_CHECK(spectrum.Evaluate(1e100) > 0.0);
    BOOST_CHECK(std::isfinite(spectrum.Evaluate(1e300)));
    auto zero = HITSpectrum::VonKarman(length, 0.0);
    BOOST_CHECK(zero.IsIdenticallyZero());
    BOOST_CHECK_EQUAL(zero.Evaluate(1.0), 0.0);
    BOOST_CHECK_THROW(HITSpectrum::VonKarman(0.0, amplitude), std::invalid_argument);
    BOOST_CHECK_THROW(HITSpectrum::VonKarman(length, -1.0), std::invalid_argument);
    BOOST_CHECK_THROW(HITSpectrum::VonKarman(std::numeric_limits<NekDouble>::infinity(),
                                           amplitude), std::invalid_argument);
    BOOST_CHECK_THROW(spectrum.Evaluate(-1.0), std::invalid_argument);
    BOOST_CHECK_THROW(spectrum.Evaluate(std::numeric_limits<NekDouble>::quiet_NaN()),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(TabulatedInterpolationAndExplicitExtensions)
{
    auto linear = HITSpectrum::Tabulated({1.0, 3.0}, {4.0, 2.0},
                                         HITInterpolation::Linear);
    BOOST_CHECK_EQUAL(linear.Evaluate(1.0), 4.0);
    BOOST_CHECK_EQUAL(linear.Evaluate(3.0), 2.0);
    Near(linear.Evaluate(2.0), 3.0);
    Near(linear.Evaluate(0.5), 0.25);
    Near(linear.Evaluate(std::sqrt(17.0)), 1.0);
    auto log = HITSpectrum::Tabulated({1.0, 4.0}, {16.0, 1.0});
    Near(log.Evaluate(2.0), 4.0);
    auto cut = HITSpectrum::Tabulated({1.0, 4.0}, {16.0, 1.0},
                                      HITInterpolation::LogLinear,
                                      HITLowWavenumber::Zero,
                                      HITHighWavenumber::Zero);
    BOOST_CHECK_EQUAL(cut.Evaluate(0.5), 0.0);
    BOOST_CHECK_EQUAL(cut.Evaluate(5.0), 0.0);
    BOOST_CHECK_EQUAL(cut.Evaluate(1.0), 16.0);
    BOOST_CHECK_EQUAL(cut.Evaluate(4.0), 1.0);
    auto zero = HITSpectrum::Tabulated({1.0, 2.0}, {0.0, 0.0},
                                      HITInterpolation::Linear,
                                      HITLowWavenumber::Zero,
                                      HITHighWavenumber::Zero);
    BOOST_CHECK(zero.IsIdenticallyZero());
    BOOST_CHECK_EQUAL(zero.Evaluate(1.5), 0.0);
}

BOOST_AUTO_TEST_CASE(GaussianTailHandlesExtremeWavenumberScales)
{
    // The dimensionless Gaussian decay is identical at either scale, even
    // though directly forming k^2 or beta in double would overflow/underflow.
    const NekDouble expected = 2.0 * std::pow(2.0, -5.0 / 3.0);
    for (NekDouble scale : {1e-200, 1e200})
    {
        const auto spectrum = HITSpectrum::Tabulated(
            {scale, 2.0 * scale}, {4.0, 2.0});
        BOOST_CHECK_CLOSE_FRACTION(spectrum.Evaluate(3.0 * scale), expected, 5e-13);
    }
}

BOOST_AUTO_TEST_CASE(InvalidTablesAndEnumValuesAreRejected)
{
    BOOST_CHECK_THROW(HITSpectrum::Tabulated({1.0}, {1.0}), std::invalid_argument);
    BOOST_CHECK_THROW(HITSpectrum::Tabulated({1.0, 2.0}, {1.0}), std::invalid_argument);
    BOOST_CHECK_THROW(HITSpectrum::Tabulated({1.0, 1.0}, {2.0, 1.0}),
                      std::invalid_argument);
    BOOST_CHECK_THROW(HITSpectrum::Tabulated({2.0, 1.0}, {2.0, 1.0}),
                      std::invalid_argument);
    BOOST_CHECK_THROW(HITSpectrum::Tabulated({0.0, 1.0}, {2.0, 1.0}),
                      std::invalid_argument);
    BOOST_CHECK_THROW(HITSpectrum::Tabulated({1.0, 2.0}, {2.0, -1.0}),
                      std::invalid_argument);
    BOOST_CHECK_THROW(HITSpectrum::Tabulated({1.0, 2.0}, {0.0, 1.0}),
                      std::invalid_argument);
    BOOST_CHECK_THROW(HITSpectrum::Tabulated({1.0, 2.0}, {1.0, 2.0}),
                      std::invalid_argument);
    BOOST_CHECK_THROW(HITSpectrum::Tabulated({1.0, 2.0}, {2.0, 1.0},
                                            static_cast<HITInterpolation>(99)),
                      std::invalid_argument);
    BOOST_CHECK_THROW(HITSpectrum::Tabulated({1.0, 2.0}, {2.0, 1.0},
                                            HITInterpolation::Linear,
                                            static_cast<HITLowWavenumber>(99)),
                      std::invalid_argument);
    BOOST_CHECK_THROW(HITSpectrum::Tabulated({1.0, 2.0}, {2.0, 1.0},
                                            HITInterpolation::Linear,
                                            HITLowWavenumber::Zero,
                                            static_cast<HITHighWavenumber>(99)),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(StrictSpectrumKindUnitsAndNumericParsing)
{
    std::istringstream valid("\n" + header + "# provenance\n1 4\n3 2\n");
    auto spectrum = HITSpectrum::ReadTable(valid, HITInterpolation::Linear);
    BOOST_CHECK_EQUAL(spectrum.GetWavenumbers().size(), 2U);
    BOOST_CHECK_EQUAL(spectrum.GetEnergies()[0], 4.0);
    for (const std::string text : {
             "1 4\n3 2\n",
             "# NEKTAR_HIT_SPECTRUM_V1 KIND=Longitudinal1D K_UNITS=rad/m E_UNITS=m3/s2\n1 4\n3 2\n",
             "# NEKTAR_HIT_SPECTRUM_V1 KIND=Shell3D K_UNITS=cm-1 E_UNITS=m3/s2\n1 4\n3 2\n",
             "# NEKTAR_HIT_SPECTRUM_V1 KIND=Shell3D K_UNITS=rad/m E_UNITS=cm3/s2\n1 4\n3 2\n"})
    {
        std::istringstream input(text);
        BOOST_CHECK_THROW(HITSpectrum::ReadTable(input), std::invalid_argument);
    }
    for (const std::string rows : {"1 4 9\n3 2\n", "1junk 4\n3 2\n",
                                   "1 nan\n3 2\n", "1 inf\n3 2\n",
                                   "1 4\n1 2\n", "1 4\n3 -2\n",
                                   "1 4\n", "1 4\n3 2 # tail\n"})
    {
        std::istringstream input(header + rows);
        BOOST_CHECK_THROW(HITSpectrum::ReadTable(input), std::invalid_argument);
    }
    BOOST_CHECK_THROW(HITSpectrum::FromFile("__nonexistent_hit_spectrum__.dat"),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(ExperimentalSpectrumRetainsMeasuredValuesAndMissingPoints)
{
    const std::string dir = HIT_TEST_FIXTURE_DIR;
    auto spectrum = HITSpectrum::FromFile(dir + "/cbc1971_theta42_spectrum_3d.dat");
    BOOST_CHECK_EQUAL(spectrum.GetWavenumbers().size(), 19U);
    BOOST_CHECK_EQUAL(spectrum.GetWavenumbers().front(), 20.0);
    BOOST_CHECK_EQUAL(spectrum.GetWavenumbers().back(), 2000.0);
    BOOST_CHECK_EQUAL(spectrum.Evaluate(50.0), 4.57e-4);
    auto downstream = HITSpectrum::FromFile(dir + "/cbc1971_theta98_spectrum_3d.dat");
    BOOST_CHECK_EQUAL(downstream.Evaluate(25.0), 1.96e-4);
    BOOST_CHECK_EQUAL(downstream.Evaluate(30.0), 1.95e-4);
    BOOST_CHECK_EQUAL(downstream.Evaluate(40.0), 2.02e-4);
    auto last = HITSpectrum::FromFile(dir + "/cbc1971_theta171_spectrum_3d.dat");
    BOOST_CHECK_EQUAL(last.GetWavenumbers().size(), 18U);
    BOOST_CHECK_EQUAL(last.GetWavenumbers().back(), 1500.0);
    BOOST_CHECK_THROW(HITSpectrum::FromFile(dir + "/cbc1971_theta42_longitudinal_1d.dat"),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(TransverseHalfLatticeAndTensorEnergyDiagnostics)
{
    HITGeneratorConfig config;
    config.lengths = {2.0 * pi, 3.0 * pi, 4.0 * pi};
    config.seed = UINT64_C(20261028);
    config.kMax = 3.0;
    config.taperStart = 2.0;
    config.counts = {10, 12, 14};
    auto spectrum = HITSpectrum::VonKarman(0.8, 2.0);
    HITGenerator generator(config, spectrum);
    const NekDouble viscosity = 1.5e-5;
    const auto diagnostic = generator.Diagnostics(viscosity);
    const auto reference = TensorMoments(config, spectrum, viscosity);
    BOOST_REQUIRE(!generator.GetModes().empty());
    BOOST_CHECK_EQUAL(generator.GetModes().size() * 2, reference.fullModes);
    BOOST_CHECK_EQUAL(diagnostic.independentModes, generator.GetModes().size());
    BOOST_CHECK_EQUAL(diagnostic.modeBytes,
                      generator.GetModes().size() * sizeof(HITFourierMode));
    BOOST_CHECK_EQUAL(HITGenerator::EstimateModeBytes(config, spectrum),
                      diagnostic.modeBytes);
    Near(diagnostic.energyScale, 1.0);
    Near(diagnostic.expectedKineticEnergy, reference.energy);
    Near(diagnostic.expectedDissipation, reference.dissipation);
    std::set<HITIndex> indices;
    NekDouble realizedK = 0.0, realizedEpsilon = 0.0;
    HITTensor realizedCovariance{};
    for (const auto &mode : generator.GetModes())
    {
        BOOST_CHECK(indices.insert(mode.index).second);
        const auto first = std::find_if(mode.index.begin(), mode.index.end(),
                                        [](int value) { return value != 0; });
        BOOST_REQUIRE(first != mode.index.end());
        BOOST_CHECK(*first > 0);
        std::complex<NekDouble> divergence{};
        NekDouble k2 = 0.0, power = 0.0;
        for (int i = 0; i < 3; ++i)
        {
            BOOST_CHECK(std::abs(mode.index[i]) < config.counts[i] / 2);
            Near(mode.wavevector[i], 2.0 * pi * mode.index[i] / config.lengths[i]);
            divergence += mode.wavevector[i] * mode.coefficient[i];
            k2 += mode.wavevector[i] * mode.wavevector[i];
            power += std::norm(mode.coefficient[i]);
            for (int j = 0; j < 3; ++j)
            {
                realizedCovariance[i][j] +=
                    2.0 * std::real(mode.coefficient[i] * std::conj(mode.coefficient[j]));
            }
        }
        BOOST_CHECK_SMALL(std::abs(divergence), 1e-12 * std::sqrt(k2 * power));
        realizedK += power;
        realizedEpsilon += 2.0 * viscosity * k2 * power;
    }
    Near(diagnostic.realizedKineticEnergy, realizedK);
    Near(diagnostic.realizedDissipation, realizedEpsilon);
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            Near(diagnostic.expectedCovariance[i][j], reference.covariance[i][j]);
            Near(diagnostic.realizedCovariance[i][j], realizedCovariance[i][j]);
        }
    }
    BOOST_CHECK_THROW(generator.Diagnostics(-1.0), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(TargetRMSNormalizesEnsembleRatherThanEachRealization)
{
    HITGeneratorConfig config;
    config.normalization = HITNormalization::TargetRMS;
    config.targetRMS = 0.222;
    auto spectrum = HITSpectrum::VonKarman(0.8, 1.0);
    HITGenerator generator(config, spectrum);
    const auto diagnostic = generator.Diagnostics();
    const auto reference = TensorMoments(config, spectrum, 0.0);
    const NekDouble targetK = 1.5 * config.targetRMS * config.targetRMS;
    Near(diagnostic.expectedKineticEnergy, targetK);
    Near(diagnostic.energyScale, targetK / reference.energy);
    NekDouble trace = 0.0;
    for (int i = 0; i < 3; ++i)
    {
        Near(diagnostic.expectedCovariance[i][i], config.targetRMS * config.targetRMS);
        trace += diagnostic.expectedCovariance[i][i];
    }
    Near(trace / 2.0, targetK);
    config.seed = UINT64_C(0x9e3779b97f4a7c15);
    const auto other = HITGenerator(config, spectrum).Diagnostics();
    BOOST_CHECK_EQUAL(other.energyScale, diagnostic.energyScale);
    BOOST_CHECK(other.realizedKineticEnergy != diagnostic.realizedKineticEnergy);
}

BOOST_AUTO_TEST_CASE(MeanOnlyHasNoModesOrDerivatives)
{
    HITGeneratorConfig config;
    config.normalization = HITNormalization::TargetRMS;
    config.targetRMS = 0.0;
    config.maxModeBytes = 0;
    auto spectrum = HITSpectrum::VonKarman(0.8, 1.0);
    HITGenerator targetZero(config, spectrum);
    BOOST_CHECK(targetZero.GetModes().empty());
    BOOST_CHECK_EQUAL(HITGenerator::EstimateModeBytes(config, spectrum), 0U);
    const auto sample = targetZero.Evaluate({0.1, 0.2, 0.3}, 0.2);
    BOOST_CHECK_EQUAL(sample.velocity[0], config.meanSpeed);
    for (int i = 0; i < 3; ++i)
    {
        if (i != 0)
        {
            BOOST_CHECK_EQUAL(sample.velocity[i], 0.0);
        }
        BOOST_CHECK_EQUAL(sample.timeDerivative[i], 0.0);
        for (int j = 0; j < 3; ++j)
        {
            BOOST_CHECK_EQUAL(sample.gradient[i][j], 0.0);
        }
    }
    BOOST_CHECK_EQUAL(targetZero.Diagnostics().expectedKineticEnergy, 0.0);
    config.normalization = HITNormalization::Preserve;
    HITGenerator preserveZero(config, HITSpectrum::VonKarman(0.8, 0.0));
    BOOST_CHECK(preserveZero.GetModes().empty());
    CheckSameSample(sample, preserveZero.Evaluate({0.1, 0.2, 0.3}, 0.2));
}

BOOST_AUTO_TEST_CASE(RealSynthesisMatchesExplicitConjugateSum)
{
    HITGeneratorConfig config;
    config.origin = {0.7, -0.4, 1.2};
    config.timeOrigin = 2.0;
    config.seed = 1234;
    HITGenerator generator(config, HITSpectrum::VonKarman(0.8, 1.0));
    const HITVector point{1.1, 0.8, -0.1};
    const NekDouble time = 2.27;
    HITComplexVector sum{};
    for (const auto &mode : generator.GetModes())
    {
        NekDouble phase = -mode.wavevector[0] * config.meanSpeed * (time - config.timeOrigin);
        for (int j = 0; j < 3; ++j)
        {
            phase += mode.wavevector[j] * (point[j] - config.origin[j]);
        }
        const auto oscillation = std::exp(std::complex<NekDouble>(0.0, phase));
        for (int i = 0; i < 3; ++i)
        {
            sum[i] += mode.coefficient[i] * oscillation +
                      std::conj(mode.coefficient[i]) * std::conj(oscillation);
        }
    }
    const auto sample = generator.Evaluate(point, time);
    for (int i = 0; i < 3; ++i)
    {
        BOOST_CHECK_SMALL(std::imag(sum[i]), 1e-14);
        Near(sample.velocity[i], std::real(sum[i]) + (i == 0 ? config.meanSpeed : 0.0));
    }
}

BOOST_AUTO_TEST_CASE(AnalyticDerivativesDivergenceAndFrozenConvection)
{
    HITGeneratorConfig config;
    config.seed = 87;
    config.meanSpeed = 1.3;
    HITGenerator generator(config, HITSpectrum::VonKarman(0.8, 1.0));
    const HITVector point{0.21, -0.43, 0.67};
    const NekDouble time = 0.2, h = 1e-5;
    const auto sample = generator.Evaluate(point, time);
    NekDouble trace = 0.0, gradientScale = 0.0;
    for (int j = 0; j < 3; ++j)
    {
        HITVector plus = point, minus = point;
        plus[j] += h;
        minus[j] -= h;
        const auto forward = generator.Evaluate(plus, time);
        const auto backward = generator.Evaluate(minus, time);
        for (int i = 0; i < 3; ++i)
        {
            Near((forward.velocity[i] - backward.velocity[i]) / (2.0 * h),
                 sample.gradient[i][j], 5e-9);
            gradientScale += std::abs(sample.gradient[i][j]);
        }
        trace += sample.gradient[j][j];
    }
    BOOST_CHECK_SMALL(trace, 1e-12 * std::max(NekDouble(1.0), gradientScale));
    const auto later = generator.Evaluate(point, time + h);
    const auto earlier = generator.Evaluate(point, time - h);
    for (int i = 0; i < 3; ++i)
    {
        Near((later.velocity[i] - earlier.velocity[i]) / (2.0 * h),
             sample.timeDerivative[i], 5e-9);
        Near(sample.timeDerivative[i], -config.meanSpeed * sample.gradient[i][0]);
    }
    HITVector convected = point;
    convected[0] += config.meanSpeed * 0.1;
    const auto followed = generator.Evaluate(convected, time + 0.1);
    for (int i = 0; i < 3; ++i)
    {
        Near(followed.velocity[i], sample.velocity[i]);
    }
}

BOOST_AUTO_TEST_CASE(AbsoluteTimeSamplingReconstructionAndPeriodicity)
{
    HITGeneratorConfig config;
    config.seed = UINT64_MAX;
    auto spectrum = HITSpectrum::VonKarman(0.8, 1.0);
    HITGenerator first(config, spectrum), reconstructed(config, spectrum);
    const HITVector point{0.2, 0.3, 0.4};
    const auto initial = first.Evaluate(point, 0.1);
    first.Evaluate({3.0, 2.0, 1.0}, 0.4);
    first.Evaluate(point, 0.0);
    CheckSameSample(initial, first.Evaluate(point, 0.1));
    CheckSameSample(initial, reconstructed.Evaluate(point, 0.1));
    for (int direction = 0; direction < 3; ++direction)
    {
        HITVector shifted = point;
        shifted[direction] += config.lengths[direction];
        const auto periodic = first.Evaluate(shifted, 0.1);
        for (int i = 0; i < 3; ++i)
        {
            Near(periodic.velocity[i], initial.velocity[i], 8e-13);
        }
    }
    config.allowRepeat = true;
    const NekDouble period = config.lengths[0] / config.meanSpeed;
    config.recordDuration = period + 0.2;
    HITGenerator repeating(config, spectrum);
    const auto afterPeriod = repeating.Evaluate(point, 0.1 + period);
    for (int i = 0; i < 3; ++i)
    {
        Near(afterPeriod.velocity[i], initial.velocity[i], 8e-13);
    }
    BOOST_CHECK_THROW(first.Evaluate(point, -0.1), std::invalid_argument);
    BOOST_CHECK_THROW(first.Evaluate(point, 0.7), std::invalid_argument);
    BOOST_CHECK_THROW(first.Evaluate({std::numeric_limits<NekDouble>::infinity(), 0.0, 0.0},
                                    0.1), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(ConfigurationAndMemoryGuards)
{
    const auto spectrum = HITSpectrum::VonKarman(0.8, 1.0);
    const HITGeneratorConfig valid;
    for (int fault = 0; fault < 9; ++fault)
    {
        auto config = valid;
        switch (fault)
        {
            case 0: config.lengths[0] = 0.0; break;
            case 1: config.counts[0] = 7; break;
            case 2: config.meanSpeed = -1.0; break;
            case 3: config.taperStart = config.kMax; break;
            case 4: config.kMax = 4.0; break;
            case 5: config.timeStep = 2.0; break;
            case 6: config.recordDuration = 7.0; break;
            case 7: config.normalization = static_cast<HITNormalization>(99); break;
            case 8: config.origin[1] = std::numeric_limits<NekDouble>::quiet_NaN(); break;
        }
        BOOST_CHECK_THROW(HITGenerator(config, spectrum), std::invalid_argument);
    }
    auto memory = valid;
    const auto bytes = HITGenerator::EstimateModeBytes(memory, spectrum);
    BOOST_REQUIRE(bytes > 0);
    memory.maxModeBytes = bytes - 1;
    BOOST_CHECK_THROW(HITGenerator(memory, spectrum), std::length_error);
    auto huge = valid;
    huge.counts = {2147483646, 2147483646, 2147483646};
    BOOST_CHECK_THROW(HITGenerator(huge, spectrum), std::overflow_error);
    auto empty = valid;
    empty.kMax = 0.9;
    empty.taperStart = 0.8;
    empty.normalization = HITNormalization::TargetRMS;
    empty.targetRMS = 0.1;
    BOOST_CHECK_THROW(HITGenerator(empty, spectrum), std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(LargeTimeOriginDoesNotWidenTheDeclaredRecord)
{
    HITGeneratorConfig config;
    config.timeOrigin = 1e16;
    config.recordDuration = 2.0;
    config.allowRepeat = true;
    HITGenerator generator(config, HITSpectrum::VonKarman(0.8, 1.0));
    const HITVector point{0.2, 0.3, 0.4};
    BOOST_CHECK_NO_THROW(generator.Evaluate(point, config.timeOrigin));
    BOOST_CHECK_NO_THROW(generator.Evaluate(point, config.timeOrigin + 2.0));
    // A tolerance scaled only by absolute time would incorrectly admit this
    // sample, fifty record lengths beyond its phase origin.
    BOOST_CHECK_THROW(generator.Evaluate(point, config.timeOrigin + 100.0),
                      std::invalid_argument);
    BOOST_CHECK_THROW(generator.Evaluate(point, config.timeOrigin - 100.0),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(UnrepresentableNormalizationScaleIsRejected)
{
    HITGeneratorConfig config;
    config.normalization = HITNormalization::TargetRMS;
    config.targetRMS = 1e-40;
    const auto excessive = HITSpectrum::Tabulated(
        {1.0, 4.0}, {1e300, 1e300}, HITInterpolation::LogLinear,
        HITLowWavenumber::K4, HITHighWavenumber::Zero);
    // The desired physical energy is representable, but the required energy
    // normalization is smaller than a double. Reporting energyScale=0 while
    // using a nonzero long-double coefficient scale would violate the contract.
    BOOST_CHECK_THROW(HITGenerator(config, excessive), std::overflow_error);

    const auto representable = HITSpectrum::Tabulated(
        {1.0, 4.0}, {1e200, 1e200}, HITInterpolation::LogLinear,
        HITLowWavenumber::K4, HITHighWavenumber::Zero);
    const auto diagnostic = HITGenerator(config, representable).Diagnostics();
    BOOST_CHECK(diagnostic.energyScale > 0.0);
    BOOST_CHECK_CLOSE_FRACTION(diagnostic.expectedKineticEnergy,
                               1.5 * config.targetRMS * config.targetRMS, 2e-13);
}

BOOST_AUTO_TEST_CASE(GaussianCoefficientEnsembleHasProperComplexCovariance)
{
    constexpr std::size_t samples = 4096;
    HITComplexVector mean{};
    HITVector fourthMoment{};
    std::array<HITComplexVector, 3> covariance{}, pseudoCovariance{};
    for (std::size_t n = 0; n < samples; ++n)
    {
        // Widely spaced seeds avoid overlapping shifted counter windows between
        // adjacent integer seeds in the frozen seed+6*p+q contract.
        const auto seed = UINT64_C(0x9e3779b97f4a7c15) * (n + 1);
        const auto gaussian = HITGenerator::Gaussian(seed, 293);
        for (int i = 0; i < 3; ++i)
        {
            mean[i] += gaussian[i] / NekDouble(samples);
            fourthMoment[i] += std::pow(std::norm(gaussian[i]), 2) /
                               NekDouble(samples);
            for (int j = 0; j < 3; ++j)
            {
                covariance[i][j] += gaussian[i] * std::conj(gaussian[j]) /
                                     NekDouble(samples);
                pseudoCovariance[i][j] += gaussian[i] * gaussian[j] /
                                          NekDouble(samples);
            }
        }
    }
    const NekDouble inverseRootN = 1.0 / std::sqrt(NekDouble(samples));
    // Six standard errors, using exact proper-complex-Gaussian second/fourth
    // moments. These account for sampling error without assuming exact variance
    // or zero covariance in one finite sample.
    for (int i = 0; i < 3; ++i)
    {
        BOOST_CHECK_SMALL(std::real(mean[i]), 6.0 * inverseRootN / std::sqrt(2.0));
        BOOST_CHECK_SMALL(std::imag(mean[i]), 6.0 * inverseRootN / std::sqrt(2.0));
        // E|xi|^4=2 and Var(|xi|^4)=20 for proper unit complex Gaussian xi.
        // This rejects fixed-magnitude random-phase substitutes that can have
        // the same second-order covariance.
        BOOST_CHECK_SMALL(fourthMoment[i] - 2.0,
                          6.0 * std::sqrt(20.0) * inverseRootN);
        for (int j = 0; j < 3; ++j)
        {
            Near(std::real(covariance[i][j]), i == j ? 1.0 : 0.0,
                 6.0 * inverseRootN);
            BOOST_CHECK_SMALL(std::imag(covariance[i][j]), 6.0 * inverseRootN);
            BOOST_CHECK_SMALL(std::real(pseudoCovariance[i][j]), 6.0 * inverseRootN);
            BOOST_CHECK_SMALL(std::imag(pseudoCovariance[i][j]), 6.0 * inverseRootN);
        }
    }
    BOOST_TEST_MESSAGE("Gaussian ensemble: " << samples << " spaced seeds; E[|xi0|^2]="
                                             << std::real(covariance[0][0]));
}

BOOST_AUTO_TEST_CASE(TransverseFieldEnsembleAndEnergyFluctuations)
{
    constexpr std::size_t samples = 4096;
    HITGeneratorConfig config;
    config.counts = {4, 4, 4};
    config.kMax = 1.9;
    config.taperStart = 1.4;
    config.normalization = HITNormalization::TargetRMS;
    config.targetRMS = 0.222;
    auto spectrum = HITSpectrum::VonKarman(0.8, 1.0);
    const auto base = TensorMoments(config, spectrum, 0.0);
    const NekDouble targetK = 1.5 * config.targetRMS * config.targetRMS;
    const auto theory = TensorMoments(config, spectrum, 0.0, targetK / base.energy);
    HITVector mean{};
    HITTensor covariance{};
    NekDouble sumK = 0.0, sumK2 = 0.0;
    for (std::size_t n = 0; n < samples; ++n)
    {
        config.seed = UINT64_C(0xd1b54a32d192ed03) * (n + 1);
        HITGenerator generator(config, spectrum);
        const auto sample = generator.Evaluate({0.21, 0.43, -0.31}, 0.2);
        HITVector fluctuation = sample.velocity;
        fluctuation[0] -= config.meanSpeed;
        for (int i = 0; i < 3; ++i)
        {
            mean[i] += fluctuation[i] / NekDouble(samples);
            for (int j = 0; j < 3; ++j)
            {
                covariance[i][j] += fluctuation[i] * fluctuation[j] /
                                     NekDouble(samples);
            }
        }
        const auto realized = generator.Diagnostics().realizedKineticEnergy;
        sumK += realized;
        sumK2 += realized * realized;
    }
    for (int i = 0; i < 3; ++i)
    {
        BOOST_CHECK_SMALL(mean[i], 6.0 * std::sqrt(theory.covariance[i][i] / samples));
        for (int j = 0; j < 3; ++j)
        {
            const NekDouble standardError =
                std::sqrt((theory.covariance[i][i] * theory.covariance[j][j] +
                           theory.covariance[i][j] * theory.covariance[i][j]) / samples);
            BOOST_CHECK_SMALL(covariance[i][j] - theory.covariance[i][j],
                              6.0 * standardError);
        }
    }
    const NekDouble meanK = sumK / samples;
    const NekDouble varianceK = (sumK2 - sumK * sumK / samples) / (samples - 1);
    BOOST_CHECK_SMALL(meanK - targetK, 6.0 * std::sqrt(theory.energyVariance / samples));
    const NekDouble varianceStandardError =
        std::sqrt(theory.fourthCumulant / samples +
                  2.0 * theory.energyVariance * theory.energyVariance / (samples - 1));
    BOOST_CHECK_SMALL(varianceK - theory.energyVariance, 6.0 * varianceStandardError);
    BOOST_CHECK(varianceK > 0.0);
    BOOST_TEST_MESSAGE("Field ensemble: " << samples << " spaced seeds; mean K=" << meanK
                      << ", target K=" << targetK << ", Var(K)=" << varianceK
                      << ", theoretical Var(K)=" << theory.energyVariance);
}

BOOST_AUTO_TEST_SUITE_END()
} // namespace Nektar::SolverUtils::HITUnitTests
