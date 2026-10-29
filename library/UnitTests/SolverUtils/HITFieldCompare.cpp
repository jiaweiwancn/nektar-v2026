///////////////////////////////////////////////////////////////////////////////
// File: HITFieldCompare.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#include "HITFieldIO.h"
#include "HITTestFields.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

using namespace Nektar;

int main(int argc, char *argv[])
{
    if (argc < 4)
    {
        std::cerr << "Usage: HITFieldCompare [session options] session.xml first.fld second.fld\n";
        return 2;
    }
    try
    {
        // Keep native options (including Tester's --use-opt-file), while the
        // last two positional arguments are fields rather than session files.
        std::vector<char *> arguments(argv, argv + argc - 2);
        arguments.push_back(nullptr);
        auto session = LibUtilities::SessionReader::CreateInstance(
            argc - 2, arguments.data());
        auto graph = SpatialDomains::MeshGraphIO::Read(session);
        auto fields = HITTestTools::Fields(session, graph);
        auto first = HITTestTools::ImportFields(session, fields, argv[argc - 2]);
        auto second = HITTestTools::ImportFields(session, fields, argv[argc - 1]);
        bool identity = true;
        for (const std::string key : {"HITInflow.IdentitySHA256", "HITInflow.ModesSHA256"})
        {
            identity &= first.metadata.count(key) && second.metadata.count(key) &&
                        first.metadata.at(key) == second.metadata.at(key);
        }
        bool sameTime = first.metadata.count("Time") && second.metadata.count("Time");
        if (sameTime)
        {
            std::size_t endA = 0, endB = 0;
            const auto a = std::stod(first.metadata.at("Time"), &endA);
            const auto b = std::stod(second.metadata.at("Time"), &endB);
            sameTime = endA == first.metadata.at("Time").size() &&
                endB == second.metadata.at("Time").size() &&
                std::isfinite(a) && std::isfinite(b) &&
                std::abs(a - b) <= 64 * std::numeric_limits<NekDouble>::epsilon() *
                                      std::max({1.0, std::abs(a), std::abs(b)});
        }
        std::cout << "HIT_FIELD_IDENTITY identical " << identity
                  << " same_time " << sameTime << '\n';
        std::cout << std::setprecision(17) << "HIT_FIELD_COMPARE {\n";
        int nonzeroVelocityComponents = 0;
        for (int variable = 0; variable < 4; ++variable)
        {
            Array<OneD, NekDouble> difference(fields[variable]->GetNpoints());
            NekDouble maximumDifference = 0.0;
            for (int point = 0; point < difference.size(); ++point)
            {
                difference[point] = first.physical[variable][point] -
                                    second.physical[variable][point];
                ASSERTL0(std::isfinite(difference[point]), "Non-finite field difference");
                maximumDifference =
                    std::max(maximumDifference, std::abs(difference[point]));
            }
            // ExpList::Linf(inarray) uses a signed maximum. A comparison error
            // needs the maximum magnitude over the complete communicator.
            session->GetComm()->AllReduce(maximumDifference,
                                          LibUtilities::ReduceMax);
            const auto firstNorm = fields[variable]->L2(first.physical[variable]);
            const auto secondNorm = fields[variable]->L2(second.physical[variable]);
            ASSERTL0(std::isfinite(firstNorm) && std::isfinite(secondNorm),
                     "Non-finite imported field norm");
            nonzeroVelocityComponents += variable < 3 && firstNorm > 0.0 && secondNorm > 0.0;
            std::cout << '"' << session->GetVariable(variable) << "\":{\"l2\":"
                      << fields[variable]->L2(difference)
                      << ",\"linf\":" << maximumDifference
                      << ",\"first_l2\":" << firstNorm
                      << ",\"second_l2\":" << secondNorm
                      << '}' << (variable == 3 ? "\n" : ",\n");
        }
        std::cout << "}\n";
        std::cout << "HIT_FIELD_NONZERO velocity_components "
                  << nonzeroVelocityComponents << '\n';
        session->Finalise();
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << std::endl;
        return 1;
    }
}
