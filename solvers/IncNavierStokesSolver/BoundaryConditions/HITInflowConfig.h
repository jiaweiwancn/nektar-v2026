///////////////////////////////////////////////////////////////////////////////
// File: HITInflowConfig.h
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_INCNAVIERSTOKES_HITINFLOWCONFIG_H
#define NEKTAR_INCNAVIERSTOKES_HITINFLOWCONFIG_H

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <SolverUtils/Turbulence/HITGenerator.h>
#include <memory>
#include <string>

namespace Nektar
{
struct HITInflowConfiguration
{
    int regionId = -1;
    SolverUtils::HITGeneratorConfig generator;
    std::shared_ptr<const SolverUtils::HITSpectrum> spectrum;
    std::string spectrumDigest;
    /// Canonical resolved settings for diagnostics and later restart metadata.
    std::string description;
};

class HITInflowConfig
{
public:
    /// Called on ALL session ranks before constructing local boundary objects.
    /// Returns null when neither HIT boundary tags nor a HITINFLOW block exists.
    /// Rank zero alone reads table bytes; all ranks agree on errors collectively.
    static std::shared_ptr<const HITInflowConfiguration> Load(
        const LibUtilities::SessionReaderSharedPtr &session);
    static std::string SHA256(const std::string &bytes);
};
} // namespace Nektar

#endif
