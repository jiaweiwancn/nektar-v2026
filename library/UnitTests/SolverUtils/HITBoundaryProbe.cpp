///////////////////////////////////////////////////////////////////////////////
// File: HITBoundaryProbe.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#include "HITTestFields.h"
#include <IncNavierStokesSolver/BoundaryConditions/IncBoundaryConditions.h>
#include <IncNavierStokesSolver/BoundaryConditions/HITInflowRestart.h>
#include <SolverUtils/Turbulence/HITPlaneEvaluator.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <utility>

using namespace Nektar;

int main(int argc, char *argv[])
{
    try
    {
        auto session = LibUtilities::SessionReader::CreateInstance(argc, argv);
        auto graph = SpatialDomains::MeshGraphIO::Read(session);
        auto fields = HITTestTools::Fields(session, graph);
        auto configuration = HITInflowConfig::Load(session);
        ASSERTL0(configuration, "HITBoundaryProbe requires HITINFLOW");
        int directLimit = -1;
        session->LoadParameter("HITProbeDirectSamples", directLimit, -1);
        ASSERTL0(directLimit >= -1, "HITProbeDirectSamples must be -1 or nonnegative");
        SolverUtils::HITGenerator reference(configuration->generator,
                                            *configuration->spectrum);
        const auto &config = configuration->generator;
        int temporalSamples = 0, outputLimit = -1;
        NekDouble sampleStart = config.timeOrigin, sampleDt = config.timeStep;
        session->LoadParameter("HITProbeTemporalSamples", temporalSamples, 0);
        session->LoadParameter("HITProbeSampleStart", sampleStart,
                               config.timeOrigin);
        session->LoadParameter("HITProbeSampleDt", sampleDt, config.timeStep);
        session->LoadParameter("HITProbeOutputPoints", outputLimit, -1);
        ASSERTL0(temporalSamples == 0 || temporalSamples >= 2,
                 "HITProbeTemporalSamples must be zero or at least two");
        ASSERTL0(outputLimit == -1 || outputLimit > 0,
                 "HITProbeOutputPoints must be -1 or positive");
        const bool temporal = temporalSamples > 0;
        auto uniformTime = [&](std::size_t index)
        {
            return static_cast<NekDouble>(
                static_cast<long double>(sampleStart) +
                static_cast<long double>(index) * sampleDt);
        };
        NekDouble sampleEnd = config.timeOrigin + config.recordDuration;
        if (temporal)
        {
            ASSERTL0(std::isfinite(sampleStart) && std::isfinite(sampleDt) &&
                     sampleDt > 0.0, "HIT probe time grid must be finite with positive dt");
            const long double end = static_cast<long double>(sampleStart) +
                static_cast<long double>(temporalSamples - 1) * sampleDt;
            ASSERTL0(std::isfinite(end) &&
                     std::abs(end) <= std::numeric_limits<NekDouble>::max(),
                     "HIT probe time-grid endpoint is not representable");
            sampleEnd = uniformTime(temporalSamples - 1);
            // Empty-point evaluation invokes reference.ValidateTime without
            // evaluating its Fourier modes or weakening the native bounds.
            const SolverUtils::HITPlaneEvaluator validator(reference, {});
            std::vector<SolverUtils::HITVector> empty;
            validator.Evaluate(sampleStart, empty);
            validator.Evaluate(sampleEnd, empty);
            ASSERTL0(uniformTime(1) > sampleStart &&
                     uniformTime(temporalSamples - 2) < sampleEnd,
                     "HIT probe sample spacing is not representable");
        }
        const auto identity = HITInflowRestart::Describe(*configuration, reference);
        IncBoundaryConditions manager;
        manager.Initialize(session, fields);
        const auto comm = session->GetComm();
        const auto stem = session->GetSessionName();
        std::ofstream csv(stem + ".probe-rank" +
                          std::to_string(comm->GetRank()) + ".csv");
        csv << std::setprecision(17)
            << "time,x,y,z,u,v,w,projected_u,projected_v,projected_w\n";
        Array<OneD, NekDouble> maximum(3, 0.0), nonzero(3, 0.0);
        NekDouble rawError = 0.0;
        long long points = 0, directPoints = 0, csvRows = 0;
        long long nativePointsPerTime = 0, exportedPointsPerTime = 0;
        std::vector<long long> exportedPerRank(comm->GetSize(), 0);
        int bad = 0, repeat = 1, repeatPerformed = 0;
        std::vector<std::vector<NekDouble>> previous;
        const std::vector<NekDouble> fractions{0.0, 0.5, 1.0, 0.5, 0.0};
        const std::size_t regularSamples = temporal ? temporalSamples : fractions.size();
        const std::size_t calls = regularSamples + (temporal ? 1 : 0);
        const std::size_t middle = temporal ? temporalSamples / 2 : 1;
        const NekDouble repeatTime = temporal ? uniformTime(middle) :
            config.timeOrigin + 0.5 * config.recordDuration;
        std::vector<int> nativeCounts(fields[0]->GetBndConditions().size(), -1);
        std::vector<std::vector<int>> selectedByBoundary(nativeCounts.size());
        std::set<SolverUtils::HITVector> uniqueExportCoordinates;
        for (std::size_t boundary = 0; boundary < nativeCounts.size(); ++boundary)
        {
            if (fields[0]->GetBndConditions()[boundary]->GetUserDefined() !=
                "HITInflow")
            {
                continue;
            }
            auto expansion = fields[0]->GetBndCondExpansions()[boundary];
            const int count = expansion->GetNpoints();
            nativeCounts[boundary] = count;
            if (outputLimit < 0)
            {
                continue;
            }
            Array<OneD, NekDouble> x(count), y(count), z(count);
            expansion->GetCoords(x, y, z);
            std::vector<std::pair<std::uint64_t, int>> candidates;
            candidates.reserve(count);
            for (int point = 0; point < count; ++point)
            {
                ASSERTL0(std::isfinite(x[point]) && std::isfinite(y[point]) &&
                         std::isfinite(z[point]), "HIT probe coordinates must be finite");
                candidates.emplace_back(SolverUtils::HITGenerator::Mix64(
                    static_cast<std::uint64_t>(point) + 0x9e3779b97f4a7c15ULL),
                    point);
            }
            std::sort(candidates.begin(), candidates.end());
            auto &selected = selectedByBoundary[boundary];
            for (const auto &candidate : candidates)
            {
                const auto point = candidate.second;
                if (uniqueExportCoordinates.insert(
                        {x[point], y[point], z[point]}).second)
                {
                    selected.push_back(point);
                    if (selected.size() == static_cast<std::size_t>(outputLimit))
                    {
                        break;
                    }
                }
            }
            // Dispersed hash selection avoids regular GLL endpoint strides.
            // Native index order makes the same selected trace stable in time.
            std::sort(selected.begin(), selected.end());
        }
        NekDouble previousTime = sampleStart;
        for (std::size_t iteration = 0; iteration < calls; ++iteration)
        {
            const bool repeatCall = temporal && iteration == regularSamples;
            const bool writeTime = !repeatCall;
            const bool saveCoefficients = iteration == middle;
            const bool checkRepeat = temporal ? repeatCall : iteration == 3;
            const NekDouble time = temporal ?
                (repeatCall ? repeatTime : uniformTime(iteration)) :
                config.timeOrigin + fractions[iteration] * config.recordDuration;
            if (temporal && iteration > 0 && !repeatCall)
            {
                ASSERTL0(time > previousTime,
                         "HIT probe sample spacing is not representable");
            }
            previousTime = time;
            manager.UpdatePrescribedVelocity(time);
            std::vector<std::vector<NekDouble>> coefficients;
            long long framePoints = 0, frameExported = 0;
            for (int boundary = 0; boundary < fields[0]->GetBndConditions().size();
                 ++boundary)
            {
                if (fields[0]->GetBndConditions()[boundary]->GetUserDefined() !=
                    "HITInflow")
                {
                    continue;
                }
                auto expansion = fields[0]->GetBndCondExpansions()[boundary];
                const int count = expansion->GetNpoints();
                ASSERTL0(count == nativeCounts[boundary],
                         "HIT probe native point set changed between calls");
                const auto &selected = selectedByBoundary[boundary];
                const int exported = outputLimit < 0 ? count :
                    static_cast<int>(selected.size());
                std::size_t nextExport = 0;
                framePoints += count;
                frameExported += exported;
                Array<OneD, NekDouble> x(count), y(count), z(count);
                Array<OneD, Array<OneD, NekDouble>> projected(3);
                expansion->GetCoords(x, y, z);
                for (int component = 0; component < 3; ++component)
                {
                    auto value = fields[component]->GetBndCondExpansions()[boundary];
                    if (saveCoefficients || checkRepeat)
                    {
                        coefficients.emplace_back(value->GetCoeffs().begin(),
                                                   value->GetCoeffs().end());
                    }
                    projected[component] = Array<OneD, NekDouble>(count);
                    const bool wave = value->GetWaveSpace();
                    value->SetWaveSpace(false);
                    value->BwdTrans(value->GetCoeffs(), projected[component]);
                    value->SetWaveSpace(wave);
                    if (iteration == 0 &&
                        fields[0]->GetExpType() == MultiRegions::e3DH1D)
                    {
                        const auto ids = value->GetZIDs();
                        const int perPlane = value->GetPlane(0)->GetNcoeffs();
                        for (int plane = 0; plane < ids.size(); ++plane)
                        {
                            if (ids[plane] < 2)
                            {
                                continue;
                            }
                            for (int coefficient = 0; coefficient < perPlane; ++coefficient)
                            {
                                const auto c = value->GetCoeffs()[plane * perPlane + coefficient];
                                nonzero[component] += c * c;
                            }
                        }
                    }
                }
                for (int point = 0; point < count; ++point)
                {
                    const bool exportPoint = writeTime && (outputLimit < 0 ||
                        (nextExport < selected.size() &&
                         point == selected[nextExport]));
                    const bool checkDirect = directLimit < 0 ||
                        (directLimit > 0 && point % std::max(1, count / directLimit) == 0 &&
                         point / std::max(1, count / directLimit) < directLimit);
                    SolverUtils::HITVector expected{};
                    if (checkDirect)
                    {
                        expected = reference.Evaluate(
                            {x[point], y[point], z[point]}, time).velocity;
                        ++directPoints;
                    }
                    if (exportPoint)
                    {
                        csv << time << ',' << x[point] << ',' << y[point] << ',' << z[point];
                    }
                    for (int component = 0; component < 3; ++component)
                    {
                        const auto raw = fields[component]->GetBndCondExpansions()[boundary]
                                             ->GetPhys()[point];
                        bad |= !std::isfinite(raw) ||
                               !std::isfinite(projected[component][point]);
                        if (checkDirect)
                        {
                            rawError = std::max(rawError,
                                std::abs(raw - expected[component]));
                        }
                        maximum[component] = std::max(maximum[component],
                            std::abs(projected[component][point] - raw));
                        if (exportPoint)
                        {
                            csv << ',' << raw;
                        }
                    }
                    if (exportPoint)
                    {
                        for (int component = 0; component < 3; ++component)
                        {
                            csv << ',' << projected[component][point];
                        }
                        csv << '\n';
                        ++csvRows;
                        ++nextExport;
                    }
                }
                points += count;
            }
            if (iteration == 0)
            {
                nativePointsPerTime = framePoints;
                exportedPointsPerTime = frameExported;
                exportedPerRank[comm->GetRank()] = frameExported;
            }
            else
            {
                ASSERTL0(framePoints == nativePointsPerTime &&
                         frameExported == exportedPointsPerTime,
                         "HIT probe point set changed between calls");
            }
            if (saveCoefficients)
            {
                previous = std::move(coefficients);
            }
            if (checkRepeat)
            {
                repeatPerformed = 1;
                if (coefficients.size() != previous.size())
                {
                    repeat = 0;
                }
                else
                {
                    for (std::size_t index = 0; index < coefficients.size(); ++index)
                    {
                        const auto &a = coefficients[index];
                        const auto &b = previous[index];
                        if (a.size() != b.size() || (!a.empty() &&
                            std::memcmp(a.data(), b.data(),
                                        a.size() * sizeof(NekDouble)) != 0))
                        {
                            repeat = 0;
                        }
                    }
                }
            }
        }
        csv.close();
        bad |= !csv;
        comm->AllReduce(bad, LibUtilities::ReduceMax);
        comm->AllReduce(points, LibUtilities::ReduceSum);
        comm->AllReduce(directPoints, LibUtilities::ReduceSum);
        comm->AllReduce(csvRows, LibUtilities::ReduceSum);
        comm->AllReduce(nativePointsPerTime, LibUtilities::ReduceSum);
        comm->AllReduce(exportedPointsPerTime, LibUtilities::ReduceSum);
        comm->AllReduce(exportedPerRank, LibUtilities::ReduceSum);
        comm->AllReduce(repeat, LibUtilities::ReduceMin);
        comm->AllReduce(repeatPerformed, LibUtilities::ReduceMin);
        comm->AllReduce(rawError, LibUtilities::ReduceMax);
        comm->AllReduce(maximum, LibUtilities::ReduceMax);
        comm->AllReduce(nonzero, LibUtilities::ReduceSum);
        ASSERTL0(!bad && repeat && repeatPerformed && points > 0,
                 "HIT boundary probe failed");
        if (comm->GetRank() == 0)
        {
            std::ofstream output(stem + ".probe.json");
            output << std::setprecision(17)
                   << "{\n\"ranks\":" << comm->GetSize()
                   << ",\n\"points\":" << points
                   << ",\n\"direct_reference_samples\":" << directPoints
                   << ",\n\"direct_reference_limit_per_boundary_rank_time\":" << directLimit
                   << ",\n\"raw_error_evaluated\":"
                   << (directPoints > 0 ? "true" : "false")
                   << ",\n\"temporal_samples\":" << temporalSamples
                   << ",\n\"time_grid\":{\"mode\":\""
                   << (temporal ? "uniform" : "legacy_fractions")
                   << "\",\"regular_samples\":" << regularSamples
                   << ",\"evaluation_calls\":" << calls
                   << ",\"start_s\":" << (temporal ? sampleStart : config.timeOrigin)
                   << ",\"sample_dt_s\":";
            if (temporal)
            {
                output << sampleDt;
            }
            else
            {
                output << "null";
            }
            output << ",\"end_s\":" << sampleEnd
                   << ",\"repeat_time_s\":" << repeatTime
                   << ",\"csv_excludes_repeat_check\":"
                   << (temporal ? "true" : "false") << '}'
                   << ",\n\"native_points_per_time\":" << nativePointsPerTime
                   << ",\n\"csv_export_limit_per_boundary_rank_time\":" << outputLimit
                   << ",\n\"csv_exported_points_per_time\":" << exportedPointsPerTime
                   << ",\n\"csv_exported_points_per_rank\":[";
            for (std::size_t rank = 0; rank < exportedPerRank.size(); ++rank)
            {
                output << (rank ? "," : "") << exportedPerRank[rank];
            }
            output << ']'
                   << ",\n\"csv_exported_rows\":" << csvRows
                   << ",\n\"csv_point_selection\":\""
                   << (outputLimit < 0 ? "all native points" :
                       "fixed lowest Mix64 native-index scores per boundary with unique coordinate triples per rank")
                   << "\""
                   << ",\n\"projection_error_reference\":\"unprojected native boundary physical values at all sampled points\""
                   << ",\n\"projection_error_scope\":\"all native inlet points at every evaluation call including repeat; not restricted to CSV selection\""
                   << ",\n\"csv_scope\":\"sampled native point traces; not an exact plane-average statistic\""
                   << ",\n\"raw_max_error\":" << rawError
                   << ",\n\"projection_max_error\":[" << maximum[0] << ','
                   << maximum[1] << ',' << maximum[2]
                   << "],\n\"nonzero_z_coeff_energy\":[" << nonzero[0] << ','
                   << nonzero[1] << ',' << nonzero[2]
                   << "],\n\"repeat_check_performed\":true"
                   << ",\n\"repeat_identical\":true,\n\"identity_sha256\":\""
                   << identity.at("HITInflow.IdentitySHA256")
                   << "\",\n\"modes_sha256\":\""
                   << identity.at("HITInflow.ModesSHA256")
                   << "\",\n\"generator_modes\":" << reference.GetModes().size() << "\n}\n";
            ASSERTL0(output.good(), "Cannot write HIT probe summary");
            std::cout << "HIT boundary probe: " << points << " samples, maximum raw error "
                      << rawError << std::endl;
            std::cout << std::setprecision(17)
                      << "HIT_PROBE raw_error " << rawError
                      << " projection_error "
                      << *std::max_element(maximum.begin(), maximum.end())
                      << " nonzero_z "
                      << (*std::min_element(nonzero.begin(), nonzero.end()) > 1e-12)
                      << " repeat " << repeat << std::endl;
        }
        session->Finalise();
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << std::endl;
        return 1;
    }
}
