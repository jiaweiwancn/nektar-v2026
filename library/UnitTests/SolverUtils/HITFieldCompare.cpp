///////////////////////////////////////////////////////////////////////////////
// File: HITFieldCompare.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#include "HITTestFields.h"
#include <LibUtilities/BasicUtils/FieldIO.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <set>

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
        ASSERTL0(session->GetComm()->GetSize() == 1,
                 "Run HITFieldCompare with one process; it imports all field partitions");
        auto fields = HITTestTools::Fields(session, graph);
        Array<OneD, Array<OneD, NekDouble>> first(4), second(4);
        for (int file = 0; file < 2; ++file)
        {
            std::vector<LibUtilities::FieldDefinitionsSharedPtr> definitions;
            std::vector<std::vector<NekDouble>> data;
            auto reader = LibUtilities::FieldIO::CreateForFile(session, argv[file + 2]);
            reader->Import(argv[file + 2], definitions, data);
            ASSERTL0(!definitions.empty() && definitions.size() == data.size(),
                     "No field data to compare");
            for (int variable = 0; variable < 4; ++variable)
            {
                std::string name = session->GetVariable(variable);
                auto field = fields[variable];
                const bool homogeneous = field->GetExpType() == MultiRegions::e3DH1D;
                const auto expectedDefinition = field->GetFieldDefinitions().front();
                std::set<std::pair<int, int>> expected, coverage;
                std::vector<int> planeIds{0};
                if (homogeneous)
                {
                    const auto ids = field->GetZIDs();
                    planeIds.assign(ids.begin(), ids.end());
                }
                for (int element = 0; element < field->GetExpSize(); ++element)
                {
                    for (int plane : planeIds)
                    {
                        expected.emplace(field->GetExp(element)->GetGeom()->GetGlobalID(), plane);
                    }
                }
                Vmath::Zero(field->GetNcoeffs(), field->UpdateCoeffs(), 1);
                bool found = false;
                for (int block = 0; block < definitions.size(); ++block)
                {
                    if (std::find(definitions[block]->m_fields.begin(),
                                  definitions[block]->m_fields.end(),
                                  name) !=
                        definitions[block]->m_fields.end())
                    {
                        found = true;
                        const auto &definition = definitions[block];
                        ASSERTL0(definition->m_numHomogeneousDir == (homogeneous ? 1 : 0) &&
                                     definition->m_basis == expectedDefinition->m_basis &&
                                     definition->m_uniOrder &&
                                     definition->m_numModes == expectedDefinition->m_numModes,
                                 "Field comparison requires matching representation and expansion modes");
                        if (homogeneous)
                        {
                            ASSERTL0(definition->m_homogeneousLengths.size() == 1 &&
                                std::abs(definition->m_homogeneousLengths[0] - field->GetHomoLen()) <
                                    1e-12 * field->GetHomoLen(), "Mismatched homogeneous period");
                        }
                        std::vector<int> storedPlanes{0};
                        if (homogeneous)
                        {
                            storedPlanes.assign(definition->m_homogeneousZIDs.begin(),
                                                 definition->m_homogeneousZIDs.end());
                        }
                        for (auto element : definition->m_elementIDs)
                        {
                            for (int plane : storedPlanes)
                            {
                                const auto key = std::make_pair(static_cast<int>(element), plane);
                                ASSERTL0(expected.count(key) && coverage.insert(key).second,
                                         "Extra or duplicate element/Fourier plane in field data");
                            }
                        }
                        field->ExtractDataToCoeffs(definitions[block], data[block],
                            name, field->UpdateCoeffs());
                    }
                }
                ASSERTL0(found, "Missing variable in comparison field");
                ASSERTL0(coverage == expected,
                         "Incomplete element/Fourier-plane coverage; import the full field directory");
                auto &values = file == 0 ? first[variable] : second[variable];
                values = Array<OneD, NekDouble>(field->GetNpoints());
                field->SetWaveSpace(false);
                field->BwdTrans(field->GetCoeffs(), values);
            }
        }
        std::cout << std::setprecision(17) << "HIT_FIELD_COMPARE {\n";
        for (int variable = 0; variable < 4; ++variable)
        {
            Array<OneD, NekDouble> difference(fields[variable]->GetNpoints());
            for (int point = 0; point < difference.size(); ++point)
            {
                ASSERTL0(std::isfinite(first[variable][point]) &&
                             std::isfinite(second[variable][point]),
                         "Non-finite comparison field");
                difference[point] = first[variable][point] - second[variable][point];
            }
            std::cout << '"' << session->GetVariable(variable) << "\":{\"l2\":"
                      << fields[variable]->L2(difference)
                      << ",\"linf\":" << fields[variable]->Linf(difference)
                      << ",\"first_l2\":" << fields[variable]->L2(first[variable])
                      << ",\"second_l2\":" << fields[variable]->L2(second[variable])
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
