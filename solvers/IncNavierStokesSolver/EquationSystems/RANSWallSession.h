///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_RANS_WALL_SESSION_H
#define NEKTAR_RANS_WALL_SESSION_H
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <SolverUtils/Turbulence/RANSWallGeometry.h>
#include <iomanip>
#include <regex>
#include <sstream>

namespace Nektar
{
inline std::string RANSWallPropertySignature(const LibUtilities::SessionReaderSharedPtr &session)
{
    std::ostringstream s; s << std::setprecision(17);
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
