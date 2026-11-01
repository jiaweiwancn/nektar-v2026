///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_RANS_WALL_SESSION_H
#define NEKTAR_RANS_WALL_SESSION_H
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <SolverUtils/Turbulence/RANSWallGeometry.h>
#include <iomanip>
#include <MultiRegions/ExpList.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <regex>
#include <sstream>
#include <cstdint>
#include <cctype>

namespace Nektar
{
inline std::string RANSStableFingerprint(const std::string &text)
{
    std::uint64_t hash = 14695981039346656037ULL;
    for (unsigned char c : text) { hash ^= c;hash *= 1099511628211ULL; }
    std::ostringstream out;out << std::hex << std::setw(16) << std::setfill('0') << hash;
    return out.str();
}

// Optional experimental inlet binding. Legacy sessions retain their original
// restart signatures. InitialConditions/file paths and NumSteps are excluded
// so full-history continuation remains possible. Inlet expressions and their
// referenced parameters are bound at runtime, not trusted from a filename.
inline std::string RANSExperimentalInletSignature(const LibUtilities::SessionReaderSharedPtr &session)
{
    if (!session->DefinesSolverInfo("RANSExperimentalSignature")) return "";
    std::ostringstream out;out << std::setprecision(17) << session->GetSolverInfo("RANSExperimentalSignature");
    for (const auto &name : session->GetVariables())
        if (name != "p")
        {
            ASSERTL0(session->DefinesFunction("RANSExperimentalInlet",name) &&
                     session->GetFunctionType("RANSExperimentalInlet",name) == LibUtilities::eFunctionTypeExpression,
                     "Experimental cases require named analytical inlet expressions");
            const auto expression = session->GetFunction("RANSExperimentalInlet",name)->GetExpression();
            out << ':' << name << '=' << expression;
            for (const auto &p : session->GetParameters())
                if (std::regex_search(expression,std::regex("\\b"+p.first+"\\b",std::regex::icase)))
                    out << ':' << p.first << '=' << p.second;
        }
    for (const auto &p : session->GetParameters())
    {
        std::string key=p.first;std::transform(key.begin(),key.end(),key.begin(),[](unsigned char c){return std::toupper(c);});
        if (key.find("RANSEXPERIMENT") == 0) out << ':' << p.first << '=' << p.second;
    }
    // Region 3 is the documented inlet of the experimental generator. Some
    // MPI partitions have no inlet region; the named expressions remain global.
    auto root=session->GetElement("NEKTAR/CONDITIONS/BOUNDARYCONDITIONS");
    auto compact=[](std::string value) {
        value.erase(std::remove_if(value.begin(),value.end(),[](unsigned char c){return std::isspace(c);}),value.end());
        return value;
    };
    for (auto region=root->FirstChildElement("REGION");region;region=region->NextSiblingElement("REGION"))
        if (std::string(region->Attribute("REF")) == "3")
            for (auto bc=region->FirstChildElement("D");bc;bc=bc->NextSiblingElement("D"))
            {
                const std::string variable=bc->Attribute("VAR");
                if (variable != "p")
                    ASSERTL0(compact(bc->Attribute("VALUE")) == compact(session->GetFunction("RANSExperimentalInlet",variable)->GetExpression()),
                             "Experimental inlet boundary differs from its bound profile");
            }
    return ":experimental-inlet-v1="+RANSStableFingerprint(out.str());
}

inline std::string RANSObservationConfiguration(const LibUtilities::SessionReaderSharedPtr &session)
{
    std::ostringstream out;out << std::setprecision(17) << "native-wall-observation-v1";
    for (const auto &name : session->GetVariables()) out << ':' << name;
    for (const auto &p : session->GetParameters())
    {
        std::string key=p.first;std::transform(key.begin(),key.end(),key.begin(),[](unsigned char c){return std::toupper(c);});
        if (key.find("RANS") == 0 || key == "KINVIS") out << ':' << p.first << '=' << p.second;
    }
    for (const std::string key : {"TurbulenceModel","RANSWallModel","RANSWallGeometry","RANSResolvedDiffusion",
        "RANSAWFWeight","RANSAWFSpatialWeight","RANSAWFFailurePolicy","AdvectionType","SpectralHPDealiasing"})
        if (session->DefinesSolverInfo(key)) out << ':' << key << '=' << session->GetSolverInfo(key);
    out << RANSExperimentalInletSignature(session);
    return "native-wall-observation-v1:"+RANSStableFingerprint(out.str());
}

inline void ShiftRANSWallPropertyQueries(const MultiRegions::ExpListSharedPtr &boundary,
    Array<OneD,NekDouble> &x,Array<OneD,NekDouble> &span)
{
    for (int e = 0; e < boundary->GetExpSize(); ++e)
    {
        const int start = boundary->GetPhys_Offset(e),count = boundary->GetExp(e)->GetTotPoints();
        for (auto *coordinate : {&x,&span})
        {
            double low = (*coordinate)[start],high = low;
            for (int q = start; q < start+count; ++q)
            { low = std::min(low,(*coordinate)[q]); high = std::max(high,(*coordinate)[q]); }
            const double tolerance = 128*std::numeric_limits<double>::epsilon()*std::max({1.,std::abs(low),std::abs(high)});
            for (int q = start; q < start+count; ++q)
                if (std::abs((*coordinate)[q]-low) <= tolerance || std::abs((*coordinate)[q]-high) <= tolerance)
                    (*coordinate)[q] += 1e-8*(.5*(low+high)-(*coordinate)[q]);
        }
    }
}

inline std::string RANSWallPropertySignature(const LibUtilities::SessionReaderSharedPtr &session)
{
    std::ostringstream s; s << std::setprecision(17);
    if (session->DefinesFunction("RANSWallProperties","roughness") || session->DefinesFunction("RANSWallProperties","z0"))
        s << ":one-sided-properties-v1";
    for (const std::string name : {"roughness","z0"})
        if (session->DefinesFunction("RANSWallProperties",name))
        {
            ASSERTL0(session->GetFunctionType("RANSWallProperties",name) == LibUtilities::eFunctionTypeExpression,
                     "RANS wall properties require stationary analytical expressions");
            const auto expression = session->GetFunction("RANSWallProperties",name)->GetExpression();
            ASSERTL0(!std::regex_search(expression,std::regex("\\bt\\b")),
                     "RANS wall properties must be independent of time");
            s << ':' << name << '=' << expression;
            // Bind parameters used by the expression, not unrelated runtime
            // controls such as NumSteps or checkpoint/output intervals.
            for (const auto &p : session->GetParameters())
                if (std::regex_search(expression,std::regex("\\b"+p.first+"\\b",std::regex::icase)))
                    s << ':' << p.first << '=' << p.second;
        }
    return s.str();
}

inline SolverUtils::RANSWallGeometry ReadRANSWallGeometry(
    const LibUtilities::SessionReaderSharedPtr &session, NekDouble origin)
{
    SolverUtils::RANSWallGeometry geometry;
    std::string kind; session->LoadSolverInfo("RANSWallGeometry",kind,"Flat");
    ASSERTL0(kind == "Flat" || kind == "Plane" || kind == "CosineHill",
             "RANSWallGeometry must be Flat, Plane or CosineHill");
    geometry.terrain = kind == "Plane" ? SolverUtils::RANSTerrain::Plane :
        (kind == "CosineHill" ? SolverUtils::RANSTerrain::CosineHill : SolverUtils::RANSTerrain::Flat);
    geometry.origin = origin;
    session->LoadParameter("RANSTerrainCentre",geometry.centre,0);
    session->LoadParameter("RANSTerrainSlope",geometry.slope,0);
    session->LoadParameter("RANSTerrainHeight",geometry.hillHeight,.04);
    session->LoadParameter("RANSTerrainHalfLength",geometry.hillHalfLength,.1);
    geometry.Validate();
    return geometry;
}

inline SolverUtils::AWFSpatialTransition ReadRANSWallTransition(
    const LibUtilities::SessionReaderSharedPtr &session,
    const SolverUtils::RANSWallGeometry &geometry)
{
    SolverUtils::AWFSpatialTransition transition;
    std::string kind; session->LoadSolverInfo("RANSAWFSpatialWeight",kind,"None");
    ASSERTL0(kind == "None" || kind == "Ramp", "RANSAWFSpatialWeight must be None or Ramp");
    transition.enabled = kind == "Ramp";
    if (transition.enabled)
    {
        session->LoadParameter("RANSAWFTransitionStart",transition.begin);
        session->LoadParameter("RANSAWFTransitionEnd",transition.end);
        transition.Validate();
        ASSERTL0(geometry.terrain != SolverUtils::RANSTerrain::CosineHill ||
                 transition.end < geometry.centre-geometry.hillHalfLength,
                 "The AWF transition must finish strictly upstream of the hill foot");
    }
    return transition;
}
}
#endif
