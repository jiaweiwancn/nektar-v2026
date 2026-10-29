///////////////////////////////////////////////////////////////////////////////
// File: HITInflowBenchmark.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/BoundaryConditions/HITInflowConfig.h>
#include <SolverUtils/Turbulence/HITPlaneEvaluator.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Nektar;
using namespace Nektar::SolverUtils;

namespace
{
using Clock = std::chrono::steady_clock;

int Count(const std::string &value)
{
    std::size_t consumed = 0;
    const auto count = std::stoi(value, &consumed);
    if (consumed != value.size() || count < 2)
    {
        throw std::invalid_argument("Plane point counts must be integers >= 2");
    }
    return count;
}

double Seconds(Clock::time_point begin, Clock::time_point end)
{
    return std::chrono::duration<double>(end - begin).count();
}

struct Moments
{
    std::size_t count = 0;
    long double mean = 0.0L, m2 = 0.0L;
    void Add(long double value)
    {
        ++count;
        const auto delta = value - mean;
        mean += delta / count;
        m2 += delta * (value - mean);
    }
    NekDouble Variance() const
    {
        return count ? static_cast<NekDouble>(m2 / count) : 0.0;
    }
};

void ArrayJSON(std::ostream &out, const std::array<NekDouble, 3> &values)
{
    out << '[' << values[0] << ',' << values[1] << ',' << values[2] << ']';
}
} // namespace

int main(int argc, char *argv[])
{
    if (argc < 4)
    {
        std::cerr << "Usage: HITInflowBenchmark [session options] session.xml "
                     "pointsPerY pointsPerZ\n";
        return 2;
    }
    try
    {
        const int ny = Count(argv[argc - 2]), nz = Count(argv[argc - 1]);
        const auto pointsCount = static_cast<std::size_t>(ny) * nz;
        ASSERTL0(pointsCount / static_cast<std::size_t>(ny) ==
                     static_cast<std::size_t>(nz), "Plane point count overflow");
        std::vector<char *> arguments(argv, argv + argc - 2);
        arguments.push_back(nullptr);
        auto session = LibUtilities::SessionReader::CreateInstance(
            argc - 2, arguments.data());
        ASSERTL0(session->GetComm()->GetSize() == 1,
                 "Run HITInflowBenchmark with one process");

        const auto configurationBegin = Clock::now();
        const auto configuration = HITInflowConfig::Load(session);
        const auto configurationEnd = Clock::now();
        ASSERTL0(configuration, "HITInflowBenchmark requires HITINFLOW");
        const auto generatorBegin = Clock::now();
        const HITGenerator generator(configuration->generator,
                                      *configuration->spectrum);
        const auto generatorEnd = Clock::now();
        const auto &config = generator.GetConfig();
        std::vector<HITVector> points;
        points.reserve(pointsCount);
        for (int y = 0; y < ny; ++y)
        {
            for (int z = 0; z < nz; ++z)
            {
                points.push_back({config.origin[0],
                    config.origin[1] + config.lengths[1] *
                        (static_cast<NekDouble>(y) / ny),
                    config.origin[2] + config.lengths[2] *
                        (static_cast<NekDouble>(z) / nz)});
            }
        }
        const auto planeBegin = Clock::now();
        const HITPlaneEvaluator plane(generator, points);
        const auto planeEnd = Clock::now();
        const NekDouble time = config.timeOrigin + 0.37 * config.recordDuration;
        constexpr int repetitions = 20;
        std::vector<HITVector> velocities;
        const auto batchBegin = Clock::now();
        for (int repeat = 0; repeat < repetitions; ++repeat)
        {
            plane.Evaluate(time, velocities);
        }
        const auto batchEnd = Clock::now();

        const std::size_t sampleCount = std::min<std::size_t>(64, points.size());
        std::vector<std::size_t> sampleIndices(sampleCount);
        for (std::size_t sample = 0; sample < sampleCount; ++sample)
        {
            // Spread across the entire plane without overflowing an index
            // multiplication for large point sets.
            sampleIndices[sample] =
                (points.size() - 1) / (sampleCount - 1) * sample +
                (points.size() - 1) % (sampleCount - 1) * sample /
                    (sampleCount - 1);
        }
        std::vector<HITVector> reference(sampleCount);
        const auto referenceBegin = Clock::now();
        for (std::size_t sample = 0; sample < sampleCount; ++sample)
        {
            reference[sample] = generator.Evaluate(
                points[sampleIndices[sample]], time).velocity;
        }
        const auto referenceEnd = Clock::now();

        std::array<NekDouble, 3> maximum{};
        std::array<Moments, 3> directMoments{}, batchMoments{};
        auto compare = [&]()
        {
            for (std::size_t sample = 0; sample < sampleCount; ++sample)
            {
                for (int component = 0; component < 3; ++component)
                {
                    const auto actual = velocities[sampleIndices[sample]][component];
                    const auto expected = reference[sample][component];
                    ASSERTL0(std::isfinite(actual) && std::isfinite(expected),
                             "Non-finite benchmark velocity");
                    maximum[component] = std::max(maximum[component],
                                                  std::abs(actual - expected));
                    // Subtract the prescribed mean before moment accumulation:
                    // u can be near 12.7 while its fluctuations are much smaller.
                    const auto mean = component == 0 ? config.meanSpeed : 0.0;
                    directMoments[component].Add(
                        static_cast<long double>(expected) - mean);
                    batchMoments[component].Add(
                        static_cast<long double>(actual) - mean);
                }
            }
        };
        compare();
        for (const auto checkTime : {config.timeOrigin,
                config.timeOrigin + 0.5 * config.recordDuration,
                config.timeOrigin + config.recordDuration})
        {
            plane.Evaluate(checkTime, velocities);
            for (std::size_t sample = 0; sample < sampleCount; ++sample)
            {
                reference[sample] = generator.Evaluate(
                    points[sampleIndices[sample]], checkTime).velocity;
            }
            compare();
        }

        constexpr NekDouble viscosity = 1.5e-5;
        const auto diagnostics = generator.Diagnostics(viscosity);
        std::array<NekDouble, 3> expectedCovariance{}, realizedCovariance{};
        std::array<NekDouble, 3> directVariance{}, batchVariance{}, varianceError{};
        for (int component = 0; component < 3; ++component)
        {
            expectedCovariance[component] = diagnostics.expectedCovariance[component][component];
            realizedCovariance[component] = diagnostics.realizedCovariance[component][component];
            directVariance[component] = directMoments[component].Variance();
            batchVariance[component] = batchMoments[component].Variance();
            varianceError[component] = std::abs(directVariance[component] -
                                               batchVariance[component]);
            for (const auto value : {expectedCovariance[component],
                    realizedCovariance[component], directVariance[component],
                    batchVariance[component], varianceError[component],
                    maximum[component]})
            {
                ASSERTL0(std::isfinite(value), "Non-finite benchmark diagnostic");
            }
        }
        std::ostringstream output;
        output.imbue(std::locale::classic());
        output << std::setprecision(17)
               << "{\"points_y\":" << ny << ",\"points_z\":" << nz
               << ",\"points\":" << points.size()
               << ",\"modes\":" << diagnostics.independentModes
               << ",\"mode_bytes\":" << diagnostics.modeBytes
               << ",\"cache_bytes\":" << plane.CacheBytes()
               << ",\"method\":\"" << plane.Method() << "\""
               << ",\"expected_kinetic_energy\":" << diagnostics.expectedKineticEnergy
               << ",\"realized_kinetic_energy\":" << diagnostics.realizedKineticEnergy
               << ",\"viscosity\":" << viscosity
               << ",\"expected_dissipation\":" << diagnostics.expectedDissipation
               << ",\"realized_dissipation\":" << diagnostics.realizedDissipation
               << ",\"expected_covariance_diagonal\":";
        ArrayJSON(output, expectedCovariance);
        output << ",\"realized_covariance_diagonal\":";
        ArrayJSON(output, realizedCovariance);
        output << ",\"configuration_seconds\":" << Seconds(configurationBegin, configurationEnd)
               << ",\"constructor_seconds\":" << Seconds(generatorBegin, generatorEnd)
               << ",\"plane_constructor_seconds\":" << Seconds(planeBegin, planeEnd)
               << ",\"batch_repetitions\":" << repetitions
               << ",\"batch_total_seconds\":" << Seconds(batchBegin, batchEnd)
               << ",\"batch_seconds\":" << Seconds(batchBegin, batchEnd) / repetitions
               << ",\"reference_seconds_per_point\":"
               << Seconds(referenceBegin, referenceEnd) / sampleCount
               << ",\"check_points\":" << sampleCount << ",\"check_times\":4"
               << ",\"max_abs_error\":" << *std::max_element(maximum.begin(), maximum.end())
               << ",\"max_abs_error_by_component\":";
        ArrayJSON(output, maximum);
        output << ",\"check_variance_reference\":";
        ArrayJSON(output, directVariance);
        output << ",\"check_variance_batch\":";
        ArrayJSON(output, batchVariance);
        output << ",\"check_variance_abs_error\":";
        ArrayJSON(output, varianceError);
        output << '}';
        std::ofstream file(session->GetSessionName() + ".benchmark.json");
        file << output.str() << '\n';
        ASSERTL0(file.good(), "Cannot write HIT inflow benchmark summary");
        std::cout << output.str() << std::endl;
        session->Finalise();
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << std::endl;
        return 1;
    }
}
