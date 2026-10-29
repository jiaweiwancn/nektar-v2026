///////////////////////////////////////////////////////////////////////////////
// File: HITFieldIO.h
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_UNITTESTS_HITFIELDIO_H
#define NEKTAR_UNITTESTS_HITFIELDIO_H

#include <LibUtilities/BasicUtils/FieldIO.h>
#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <MultiRegions/ExpList.h>

#include <algorithm>
#include <cmath>
#include <locale>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace Nektar::HITTestTools
{
struct ImportedField
{
    Array<OneD, Array<OneD, NekDouble>> physical;
    LibUtilities::FieldMetaDataMap metadata;
};

/// Import every field partition into a serial field set, checking full coverage.
/// Fields retain the imported coefficients with WaveSpace(false); physical data
/// are returned in separate arrays rather than written into field->GetPhys().
inline ImportedField ImportFields(
    const LibUtilities::SessionReaderSharedPtr &session,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const std::string &filename)
{
    ASSERTL0(session && session->GetComm()->GetSize() == 1,
             "Run HIT field import with one process; it imports all field partitions");
    ASSERTL0(fields.size() == 4 && session->GetVariables().size() == 4,
             "HIT field import requires four velocity/pressure fields");

    ImportedField result;
    result.physical = Array<OneD, Array<OneD, NekDouble>>(4);
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> definitions;
    std::vector<std::vector<NekDouble>> data;
    auto reader = LibUtilities::FieldIO::CreateForFile(session, filename);
    // Import into a fresh map: native imports merge rather than clear metadata.
    reader->Import(filename, definitions, data, result.metadata);
    ASSERTL0(!definitions.empty() && definitions.size() == data.size(),
             "No field data to compare");

    for (int variable = 0; variable < 4; ++variable)
    {
        std::string name = session->GetVariable(variable);
        auto field = fields[variable];
        ASSERTL0(field, "Missing expansion in HIT field import");
        const bool homogeneous =
            field->GetExpType() == MultiRegions::e3DH1D;
        const auto expectedDefinitions = field->GetFieldDefinitions();
        ASSERTL0(!expectedDefinitions.empty(),
                 "No expansion definitions for HIT field import");
        const auto expectedDefinition = expectedDefinitions.front();
        // GeneralGetFieldDefinitions retains one mode tuple per element even
        // for uniform order; FieldIO writes only the first UNIORDER tuple.
        const auto dimensions = expectedDefinition->m_basis.size();
        ASSERTL0(expectedDefinition->m_uniOrder && dimensions > 0 &&
                     expectedDefinition->m_numModes.size() >= dimensions,
                 "Field comparison requires uniform expansion modes");

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
                expected.emplace(
                    field->GetExp(element)->GetGeom()->GetGlobalID(), plane);
            }
        }

        Vmath::Zero(field->GetNcoeffs(), field->UpdateCoeffs(), 1);
        bool found = false;
        for (int block = 0; block < definitions.size(); ++block)
        {
            auto &definition = definitions[block];
            if (std::find(definition->m_fields.begin(),
                          definition->m_fields.end(), name) ==
                definition->m_fields.end())
            {
                continue;
            }
            found = true;
            ASSERTL0(
                definition->m_numHomogeneousDir == (homogeneous ? 1 : 0) &&
                    definition->m_basis == expectedDefinition->m_basis &&
                    definition->m_uniOrder &&
                    definition->m_numModes.size() == dimensions &&
                    std::equal(definition->m_numModes.begin(),
                               definition->m_numModes.end(),
                               expectedDefinition->m_numModes.begin()),
                "Field comparison requires matching representation and expansion modes");
            if (homogeneous)
            {
                // The native XML writer uses six significant digits for
                // HOMOGENEOUSLENGTHS. Accept that exact round trip as well as a
                // full-precision header.
                std::ostringstream xmlLength;
                xmlLength.imbue(std::locale::classic());
                xmlLength << field->GetHomoLen();
                const auto xmlPeriod = std::stod(xmlLength.str());
                ASSERTL0(
                    definition->m_homogeneousLengths.size() == 1 &&
                        (std::abs(definition->m_homogeneousLengths[0] -
                                  field->GetHomoLen()) <
                             1e-12 * field->GetHomoLen() ||
                         std::abs(definition->m_homogeneousLengths[0] -
                                  xmlPeriod) <
                             1e-12 * field->GetHomoLen()),
                    "Mismatched homogeneous period");
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
                    const auto key =
                        std::make_pair(static_cast<int>(element), plane);
                    ASSERTL0(expected.count(key) &&
                                 coverage.insert(key).second,
                             "Extra or duplicate element/Fourier plane in field data");
                }
            }
            field->ExtractDataToCoeffs(definition, data[block], name,
                                       field->UpdateCoeffs());
        }
        ASSERTL0(found, "Missing variable in comparison field");
        ASSERTL0(coverage == expected,
                 "Incomplete element/Fourier-plane coverage; import the full field directory");

        auto &values = result.physical[variable];
        values = Array<OneD, NekDouble>(field->GetNpoints());
        field->SetWaveSpace(false);
        field->BwdTrans(field->GetCoeffs(), values);
        for (const auto value : values)
        {
            ASSERTL0(std::isfinite(value), "Non-finite comparison field");
        }
    }
    return result;
}
} // namespace Nektar::HITTestTools

#endif
