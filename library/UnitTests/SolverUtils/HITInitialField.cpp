///////////////////////////////////////////////////////////////////////////////
// File: HITInitialField.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#include "HITTestFields.h"
#include <IncNavierStokesSolver/BoundaryConditions/HITInflowRestart.h>
#include <IncNavierStokesSolver/BoundaryConditions/IncBoundaryConditions.h>
#include <SolverUtils/Turbulence/HITPlaneEvaluator.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <map>
#include <memory>
#include <sstream>
#include <vector>

using namespace Nektar;

namespace
{
using Clock = std::chrono::steady_clock;

NekDouble Seconds(const Clock::time_point &start)
{
    return std::chrono::duration<NekDouble>(Clock::now() - start).count();
}

std::string Number(NekDouble value)
{
    std::ostringstream output;
    output.imbue(std::locale::classic());
    output << std::setprecision(std::numeric_limits<NekDouble>::max_digits10)
           << value;
    return output.str();
}

// Sampling has no collectives. Agree on local errors before peers enter the
// collective transforms or field writer, including ranks without inlet faces.
void Agree(const LibUtilities::CommSharedPtr &comm, std::string error)
{
    int owner = error.empty() ? comm->GetSize() : comm->GetRank();
    comm->AllReduce(owner, LibUtilities::ReduceMin);
    if (owner == comm->GetSize())
    {
        return;
    }
    int length = comm->GetRank() == owner ? error.size() : 0;
    comm->Bcast(length, owner);
    std::vector<char> bytes(length);
    if (comm->GetRank() == owner)
    {
        std::copy(error.begin(), error.end(), bytes.begin());
    }
    comm->Bcast(bytes, owner);
    NEKERROR(ErrorUtil::efatal,
             "HITInitialField: " + std::string(bytes.begin(), bytes.end()));
}
} // namespace

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        std::cerr << "Usage: HITInitialField [session options] session.xml output.fld\n";
        return 2;
    }
    try
    {
        const auto setupStart = Clock::now();
        const std::string outputFilename = argv[argc - 1];
        std::vector<char *> arguments(argv, argv + argc - 1);
        arguments.push_back(nullptr);
        auto session = LibUtilities::SessionReader::CreateInstance(
            argc - 1, arguments.data());
        auto graph = SpatialDomains::MeshGraphIO::Read(session);
        auto fields = HITTestTools::Fields(session, graph);
        const auto comm = session->GetComm();
        const auto configuration = HITInflowConfig::Load(session);
        ASSERTL0(configuration, "HITInitialField requires HITINFLOW");

        std::unique_ptr<SolverUtils::HITGenerator> generator;
        std::string error;
        try
        {
            generator = std::make_unique<SolverUtils::HITGenerator>(
                configuration->generator, *configuration->spectrum);
        }
        catch (const std::exception &exception)
        {
            error = exception.what();
        }
        Agree(comm, error);

        const NekDouble time = configuration->generator.timeOrigin;
        IncBoundaryConditions manager;
        manager.Initialize(session, fields);
        manager.UpdatePrescribedVelocity(time);
        NekDouble setupSeconds = Seconds(setupStart);

        const int count = fields[0]->GetNpoints();
        unsigned long long globalPoints = count, globalGroups = 0;
        unsigned long long grouped = 0, reference = 0, cachePeak = 0;
        unsigned long long groupsMaximum = 0;
        const auto samplingStart = Clock::now();
        error.clear();
        try
        {
            Array<OneD, NekDouble> x(count), y(count), z(count);
            fields[0]->GetCoords(x, y, z);
            std::map<NekDouble, std::vector<int>> groups;
            for (int point = 0; point < count; ++point)
            {
                ASSERTL0(std::isfinite(x[point]) && std::isfinite(y[point]) &&
                             std::isfinite(z[point]),
                         "Non-finite initial-field coordinate");
                groups[x[point]].push_back(point);
            }
            globalGroups = groupsMaximum = groups.size();
            // Keep only one plane evaluator/cache alive. The groups contain
            // indices, while the resulting velocity is written into the fields.
            for (const auto &group : groups)
            {
                std::vector<SolverUtils::HITVector> points;
                points.reserve(group.second.size());
                for (int point : group.second)
                {
                    points.push_back({x[point], y[point], z[point]});
                }
                SolverUtils::HITPlaneEvaluator evaluator(*generator, points);
                std::vector<SolverUtils::HITVector> velocity;
                evaluator.Evaluate(time, velocity);
                ASSERTL0(velocity.size() == group.second.size(),
                         "Initial-field evaluator returned the wrong sample count");
                grouped += evaluator.IsGrouped();
                reference += !evaluator.IsGrouped();
                cachePeak = std::max(cachePeak,
                    static_cast<unsigned long long>(evaluator.CacheBytes()));
                for (std::size_t sample = 0; sample < group.second.size(); ++sample)
                {
                    const int point = group.second[sample];
                    for (int component = 0; component < 3; ++component)
                    {
                        ASSERTL0(std::isfinite(velocity[sample][component]),
                                 "Non-finite initial velocity");
                        fields[component]->UpdatePhys()[point] =
                            velocity[sample][component];
                    }
                    fields[3]->UpdatePhys()[point] = 0.0;
                }
            }
        }
        catch (const std::exception &exception)
        {
            error = exception.what();
        }
        Agree(comm, error);
        NekDouble samplingSeconds = Seconds(samplingStart);

        const auto projectionStart = Clock::now();
        for (int variable = 0; variable < 4; ++variable)
        {
            auto field = fields[variable];
            field->SetWaveSpace(false);
            // For Homogeneous1D this performs the homogeneous forward
            // transform exactly once, after constrained planar projection.
            field->FwdTransBndConstrained(field->GetPhys(), field->UpdateCoeffs());
            field->AvgAssemble();
            field->ImposeDirichletConditions(field->UpdateCoeffs());
            field->LocalToGlobal();
            field->GlobalToLocal();
            field->SetPhysState(false);
        }
        NekDouble projectionSeconds = Seconds(projectionStart);

        NekDouble viscosity = 0.0;
        session->LoadParameter("Kinvis", viscosity, 0.0);
        const auto diagnostics = generator->Diagnostics(viscosity);
        auto metadata = HITInflowRestart::Describe(*configuration, *generator);
        metadata["Time"] = Number(time);
        metadata["TimeStep"] = Number(configuration->generator.timeStep);
        metadata["Kinvis"] = Number(viscosity);
        metadata["ChkFileNum"] = "0";
        // This is a generated initial condition, not an imported restart:
        // Describe supplies identity without attempting to read InitialConditions.
        auto definitions = fields[0]->GetFieldDefinitions();
        std::vector<std::vector<NekDouble>> data(definitions.size());
        for (int variable = 0; variable < 4; ++variable)
        {
            for (std::size_t block = 0; block < definitions.size(); ++block)
            {
                definitions[block]->m_fields.push_back(session->GetVariable(variable));
                fields[0]->AppendFieldData(definitions[block], data[block],
                                           fields[variable]->UpdateCoeffs());
            }
        }
        const auto writeStart = Clock::now();
        auto writer = LibUtilities::FieldIO::CreateDefault(session);
        writer->Write(outputFilename, definitions, data, metadata,
                       session->GetBackups());
        NekDouble writeSeconds = Seconds(writeStart);

        comm->AllReduce(globalPoints, LibUtilities::ReduceSum);
        comm->AllReduce(globalGroups, LibUtilities::ReduceSum);
        comm->AllReduce(grouped, LibUtilities::ReduceSum);
        comm->AllReduce(reference, LibUtilities::ReduceSum);
        comm->AllReduce(groupsMaximum, LibUtilities::ReduceMax);
        comm->AllReduce(cachePeak, LibUtilities::ReduceMax);
        for (auto value : {&setupSeconds, &samplingSeconds,
                           &projectionSeconds, &writeSeconds})
        {
            comm->AllReduce(*value, LibUtilities::ReduceMax);
        }
        if (comm->GetRank() == 0)
        {
            std::cout << std::setprecision(17)
                << "HIT_INITIAL_FIELD {\"time\":" << time
                << ",\"ranks\":" << comm->GetSize()
                << ",\"global_quadrature_points\":" << globalPoints
                << ",\"plane_groups_sum\":" << globalGroups
                << ",\"plane_groups_max_per_rank\":" << groupsMaximum
                << ",\"grouped_planes\":" << grouped
                << ",\"reference_planes\":" << reference
                << ",\"evaluator_cache_payload_peak_bytes_per_rank\":" << cachePeak
                << ",\"mode_bytes_per_generator\":" << diagnostics.modeBytes
                << ",\"independent_modes\":" << diagnostics.independentModes
                << ",\"expected_periodic_K\":" << diagnostics.expectedKineticEnergy
                << ",\"realized_periodic_K\":" << diagnostics.realizedKineticEnergy
                << ",\"expected_periodic_dissipation\":" << diagnostics.expectedDissipation
                << ",\"realized_periodic_dissipation\":" << diagnostics.realizedDissipation
                << ",\"energy_scale\":" << diagnostics.energyScale
                << ",\"setup_seconds_max\":" << setupSeconds
                << ",\"sampling_seconds_max\":" << samplingSeconds
                << ",\"projection_seconds_max\":" << projectionSeconds
                << ",\"write_seconds_max\":" << writeSeconds << "}\n";
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
