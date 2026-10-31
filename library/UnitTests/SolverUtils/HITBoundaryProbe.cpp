///////////////////////////////////////////////////////////////////////////////
// File: HITBoundaryProbe.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#include "HITTestFields.h"
#include <IncNavierStokesSolver/BoundaryConditions/IncBoundaryConditions.h>
#include <IncNavierStokesSolver/BoundaryConditions/HITInflowRestart.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

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
        int points = 0, directPoints = 0, bad = 0, repeat = 1;
        std::vector<std::vector<NekDouble>> previous;
        const auto &config = configuration->generator;
        const std::vector<NekDouble> fractions{0.0, 0.5, 1.0, 0.5, 0.0};
        for (int iteration = 0; iteration < fractions.size(); ++iteration)
        {
            const NekDouble time = config.timeOrigin +
                                   fractions[iteration] * config.recordDuration;
            manager.UpdatePrescribedVelocity(time);
            std::vector<std::vector<NekDouble>> coefficients;
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
                Array<OneD, NekDouble> x(count), y(count), z(count);
                Array<OneD, Array<OneD, NekDouble>> projected(3);
                expansion->GetCoords(x, y, z);
                for (int component = 0; component < 3; ++component)
                {
                    auto value = fields[component]->GetBndCondExpansions()[boundary];
                    coefficients.emplace_back(value->GetCoeffs().begin(),
                                               value->GetCoeffs().end());
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
                    csv << time << ',' << x[point] << ',' << y[point] << ',' << z[point];
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
                        csv << ',' << raw;
                    }
                    for (int component = 0; component < 3; ++component)
                    {
                        csv << ',' << projected[component][point];
                    }
                    csv << '\n';
                }
                points += count;
            }
            if (iteration == 1)
            {
                previous = coefficients;
            }
            if (iteration == 3 && coefficients != previous)
            {
                repeat = 0;
            }
        }
        csv.close();
        bad |= !csv;
        comm->AllReduce(bad, LibUtilities::ReduceMax);
        comm->AllReduce(points, LibUtilities::ReduceSum);
        comm->AllReduce(directPoints, LibUtilities::ReduceSum);
        comm->AllReduce(repeat, LibUtilities::ReduceMin);
        comm->AllReduce(rawError, LibUtilities::ReduceMax);
        comm->AllReduce(maximum, LibUtilities::ReduceMax);
        comm->AllReduce(nonzero, LibUtilities::ReduceSum);
        ASSERTL0(!bad && repeat && points > 0, "HIT boundary probe failed");
        if (comm->GetRank() == 0)
        {
            std::ofstream output(stem + ".probe.json");
            output << std::setprecision(17)
                   << "{\n\"ranks\":" << comm->GetSize()
                   << ",\n\"points\":" << points
                   << ",\n\"direct_reference_samples\":" << directPoints
                   << ",\n\"direct_reference_limit_per_boundary_rank_time\":" << directLimit
                   << ",\n\"projection_error_reference\":\"unprojected native boundary physical values at all sampled points\""
                   << ",\n\"raw_max_error\":" << rawError
                   << ",\n\"projection_max_error\":[" << maximum[0] << ','
                   << maximum[1] << ',' << maximum[2]
                   << "],\n\"nonzero_z_coeff_energy\":[" << nonzero[0] << ','
                   << nonzero[1] << ',' << nonzero[2]
                   << "],\n\"repeat_identical\":true,\n\"identity_sha256\":\""
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
