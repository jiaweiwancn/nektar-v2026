///////////////////////////////////////////////////////////////////////////////
// File: RANSFieldAudit.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <LibUtilities/BasicUtils/FieldIO.h>
#include <MultiRegions/ContField.h>
#include <SpatialDomains/MeshGraphIO.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <fstream>

using namespace Nektar;
using Fields = Array<OneD, MultiRegions::ExpListSharedPtr>;
using Values = Array<OneD, Array<OneD, NekDouble>>;

Values Import(const LibUtilities::SessionReaderSharedPtr &session,
              const Fields &fields, const std::vector<std::string> &names, const std::string &path,
              LibUtilities::FieldMetaDataMap &metadata)
{
    std::vector<LibUtilities::FieldDefinitionsSharedPtr> definitions;
    std::vector<std::vector<NekDouble>> data;
    LibUtilities::FieldIO::CreateForFile(session, path)->Import(path, definitions, data, metadata);
    ASSERTL0(!definitions.empty(), "Empty RANS field import");
    Values values(fields.size());
    for (int f = 0; f < fields.size(); ++f)
    {
        std::set<int> expected, seen;
        for (int e = 0; e < fields[f]->GetExpSize(); ++e)
            expected.insert(fields[f]->GetExp(e)->GetGeom()->GetGlobalID());
        Vmath::Zero(fields[f]->GetNcoeffs(), fields[f]->UpdateCoeffs(), 1);
        std::string name = names[f];
        for (int b = 0; b < definitions.size(); ++b)
        {
            auto &definition = definitions[b];
            if (std::find(definition->m_fields.begin(), definition->m_fields.end(), name) == definition->m_fields.end())
                continue;
            ASSERTL0(definition->m_numHomogeneousDir == 0, "Full-dimensional fields required");
            const auto expectedDefinition = fields[f]->GetFieldDefinitions().front();
            ASSERTL0(definition->m_uniOrder && definition->m_basis == expectedDefinition->m_basis &&
                         definition->m_numModes.size() == expectedDefinition->m_basis.size() &&
                         std::equal(definition->m_numModes.begin(), definition->m_numModes.end(),
                                    expectedDefinition->m_numModes.begin()),
                     "Mismatched imported expansion modes");
            for (int id : definition->m_elementIDs)
                ASSERTL0(expected.count(id) && seen.insert(id).second, "Duplicate or extra field element");
            fields[f]->ExtractDataToCoeffs(definition, data[b], name, fields[f]->UpdateCoeffs());
        }
        ASSERTL0(seen == expected, "Incomplete field coverage; import the complete MPI field directory");
        values[f] = Array<OneD, NekDouble>(fields[f]->GetTotPoints());
        fields[f]->BwdTrans(fields[f]->GetCoeffs(), values[f]);
    }
    return values;
}

int main(int argc, char **argv)
{
    if (argc != 4 && argc != 5)
    {
        std::cerr << "Usage: RANSFieldAudit session.xml first.fld second.fld [profiles.csv]\n";
        return 2;
    }
    try
    {
        char *arguments[]{argv[0], argv[1], nullptr};
        auto session = LibUtilities::SessionReader::CreateInstance(2, arguments);
        ASSERTL0(session->GetComm()->GetSize() == 1, "Audit all partitions with one process");
        auto graph = SpatialDomains::MeshGraphIO::Read(session);
        std::vector<std::string> names = session->GetVariables();
        for (const auto &name : {"rans_nut", "rans_P", "rans_Yap", "rans_eps_total",
             "rans_R_xx", "rans_R_yy", "rans_R_zz", "rans_R_xy", "rans_R_xz", "rans_R_yz",
             "rans_quadratic_norm", "rans_cubic_norm"})
            if (session->DefinesFunction("ExactSolution", name)) names.push_back(name);
        Fields fields(names.size());
        // EquationSystem::WriteFld embeds lower-order pressure coefficients
        // in the first field's expansion. Audit that complete saved space,
        // including the padded modes, before comparing physical values.
        fields[0] = MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(
            session, graph, session->GetVariable(0), true, false);
        for (int f = 1; f < fields.size(); ++f) fields[f] = fields[0];
        LibUtilities::FieldMetaDataMap firstMeta, secondMeta;
        auto first = Import(session, fields, names, argv[2], firstMeta);
        auto second = Import(session, fields, names, argv[3], secondMeta);
        ASSERTL0(firstMeta.count("Time") && secondMeta.count("Time"), "Missing field time");
        const double time = std::stod(firstMeta.at("Time"));
        ASSERTL0(std::abs(time - std::stod(secondMeta.at("Time"))) < 1e-12, "Mismatched field times");
        const int nq = fields[0]->GetTotPoints();
        Array<OneD, NekDouble> x(nq, 0.0), y(nq, 0.0), z(nq, 0.0), exact(nq), difference(nq);
        fields[0]->GetCoords(x, y, z);
        std::cout << std::setprecision(17) << "{\"time\":" << time << ",\"fields\":{";
        for (int f = 0; f < fields.size(); ++f)
        {
            const auto &name = names[f];
            session->GetFunction("ExactSolution", name)->Evaluate(x, y, z, time, exact);
            double errorMax = 0, differenceMax = 0;
            double minimum = std::numeric_limits<double>::max(), maximum = -minimum;
            for (int q = 0; q < nq; ++q)
            {
                ASSERTL0(std::isfinite(first[f][q]) && std::isfinite(second[f][q]), "Nonfinite field value");
                difference[q] = first[f][q] - exact[q];
                errorMax = std::max(errorMax, std::abs(difference[q]));
                differenceMax = std::max(differenceMax, std::abs(first[f][q] - second[f][q]));
                minimum = std::min(minimum, first[f][q]); maximum = std::max(maximum, first[f][q]);
            }
            const double errorL2 = fields[f]->L2(difference);
            Vmath::Vsub(nq, first[f], 1, second[f], 1, difference, 1);
            if (f) std::cout << ',';
            std::cout << '"' << name << "\":{\"l2_error\":" << errorL2
                      << ",\"linf_error\":" << errorMax << ",\"difference_l2\":" << fields[f]->L2(difference)
                      << ",\"difference_linf\":" << differenceMax << ",\"min\":" << minimum
                      << ",\"max\":" << maximum << '}';
        }
        std::cout << '}';
        if (firstMeta.count("RANSHistoryFile"))
            std::cout << ",\"history_directory\":" << std::quoted(firstMeta.at("RANSHistoryFile"));
        std::cout << "}\n";
        if (argc == 5)
        {
            std::ofstream csv(argv[4]);
            ASSERTL0(csv.good(), "Cannot open profile output");
            csv << std::setprecision(17) << "x,height,field,value,reference\n";
            double length, height, top, span;
            session->LoadParameter("L", length); session->LoadParameter("H", height);
            session->LoadParameter("Z", top); session->LoadParameter("Span", span, 1.0);
            const int dim = fields[0]->GetCoordim(0);
            std::vector<double> heights;
            for (int j = 0; j <= 64; ++j)
                heights.push_back(j == 64 ? top : height*std::pow(top/height,j/64.0));
            // Common physical samples also permit comparison when H changes.
            for (double value : {4.0,8.0,10.0,16.0,20.0,25.0,30.0})
                if (value >= height && value <= top) heights.push_back(value);
            std::sort(heights.begin(),heights.end());
            heights.erase(std::unique(heights.begin(),heights.end(),[](double a,double b) {
                return std::abs(a-b)<1e-10;
            }),heights.end());
            for (double station : {.25, .5, .75, .95})
                for (double altitude : heights)
                {
                    Array<OneD, NekDouble> point(dim), local(dim);
                    point[0] = station*length; point[1] = dim == 2 ? altitude : span*.5;
                    if (dim == 3) point[2] = altitude;
                    const int e = fields[0]->GetExpIndex(point, local, 1e-8);
                    ASSERTL0(e >= 0, "Profile point outside mesh");
                    Array<OneD, NekDouble> px(1,point[0]), py(1,point[1]), pz(1,dim == 3 ? point[2] : 0.0), value(1);
                    for (int f = 0; f < fields.size(); ++f)
                    {
                        session->GetFunction("ExactSolution",names[f])->Evaluate(px,py,pz,time,value);
                        csv << point[0] << ',' << altitude << ',' << names[f] << ','
                            // GetExpIndex already returned element-local coordinates.
                            << fields[f]->GetExp(e)->StdPhysEvaluate(local,first[f]+fields[f]->GetPhys_Offset(e))
                            << ',' << value[0] << '\n';
                    }
                }
        }
        session->Finalise();
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n'; return 1;
    }
}
