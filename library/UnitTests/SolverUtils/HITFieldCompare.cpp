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

using namespace Nektar;

int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        std::cerr << "Usage: HITFieldCompare session.xml first.fld second.fld\n";
        return 2;
    }
    try
    {
        char *arguments[] = {argv[0], argv[1], nullptr};
        auto session = LibUtilities::SessionReader::CreateInstance(2, arguments);
        auto graph = SpatialDomains::MeshGraphIO::Read(session);
        auto fields = HITTestTools::Fields(session, graph);
        auto first = HITTestTools::ImportFields(session, fields, argv[2]);
        auto second = HITTestTools::ImportFields(session, fields, argv[3]);
        bool identity = true;
        for (const std::string key : {"HITInflow.IdentitySHA256", "HITInflow.ModesSHA256"})
        {
            identity &= first.metadata.count(key) && second.metadata.count(key) &&
                        first.metadata.at(key) == second.metadata.at(key);
        }
        bool sameTime = first.metadata.count("Time") && second.metadata.count("Time");
        if (sameTime)
        {
            const auto a = std::stod(first.metadata.at("Time"));
            const auto b = std::stod(second.metadata.at("Time"));
            sameTime = std::isfinite(a) && std::isfinite(b) &&
                std::abs(a - b) <= 64 * std::numeric_limits<NekDouble>::epsilon() *
                                      std::max({1.0, std::abs(a), std::abs(b)});
        }
        std::cout << "HIT_FIELD_IDENTITY identical " << identity
                  << " same_time " << sameTime << '\n';
        std::cout << std::setprecision(17) << "HIT_FIELD_COMPARE {\n";
        for (int variable = 0; variable < 4; ++variable)
        {
            Array<OneD, NekDouble> difference(fields[variable]->GetNpoints());
            for (int point = 0; point < difference.size(); ++point)
            {
                difference[point] = first.physical[variable][point] -
                                    second.physical[variable][point];
            }
            std::cout << '"' << session->GetVariable(variable) << "\":{\"l2\":"
                      << fields[variable]->L2(difference)
                      << ",\"linf\":" << fields[variable]->Linf(difference)
                      << ",\"first_l2\":" << fields[variable]->L2(first.physical[variable])
                      << ",\"second_l2\":" << fields[variable]->L2(second.physical[variable])
                      << '}' << (variable == 3 ? "\n" : ",\n");
        }
        std::cout << "}\n";
        session->Finalise();
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << std::endl;
        return 1;
    }
}
