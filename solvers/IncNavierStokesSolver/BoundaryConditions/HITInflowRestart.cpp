///////////////////////////////////////////////////////////////////////////////
// File: HITInflowRestart.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/BoundaryConditions/HITInflowRestart.h>
#include <LibUtilities/BasicConst/GitRevision.h>
#include <algorithm>
#include <cfenv>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>

namespace Nektar
{
namespace
{
std::string Number(NekDouble value)
{
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << std::setprecision(std::numeric_limits<NekDouble>::max_digits10)
        << (value == 0.0 ? 0.0 : value);
    return out.str();
}

template <typename Vector> std::string Tuple(const Vector &values)
{
    std::ostringstream out;
    out.imbue(std::locale::classic());
    bool first = true;
    for (auto value : values)
    {
        if (!first)
        {
            out << ' ';
        }
        out << Number(value);
        first = false;
    }
    return out.str();
}

std::string Environment()
{
    std::ostringstream out;
    out.imbue(std::locale::classic());
#if defined(__VERSION__)
    out << "compiler=" << __VERSION__;
#elif defined(_MSC_VER)
    out << "compiler=MSVC-" << _MSC_VER;
#else
    out << "compiler=unknown";
#endif
#if defined(__linux__)
    out << ";platform=Linux";
#elif defined(__APPLE__)
    out << ";platform=Darwin";
#elif defined(_WIN32)
    out << ";platform=Windows";
#else
    out << ";platform=unknown";
#endif
#if defined(__GLIBC__)
    out << ";glibc=" << __GLIBC__ << '.' << __GLIBC_MINOR__;
#endif
    out << ";cplusplus=" << __cplusplus << ";sizeof=" << sizeof(NekDouble)
        << ";digits=" << std::numeric_limits<NekDouble>::digits
        << ";iec559=" << std::numeric_limits<NekDouble>::is_iec559
        << ";round=" << std::fegetround();
#if defined(__FAST_MATH__)
    out << ";fastmath=1";
#else
    out << ";fastmath=0";
#endif
    return out.str();
}

void Agree(const LibUtilities::CommSharedPtr &comm, std::string error)
{
    int rank = error.empty() ? comm->GetSize() : comm->GetRank();
    comm->AllReduce(rank, LibUtilities::ReduceMin);
    if (rank == comm->GetSize())
    {
        return;
    }
    int count = error.size();
    comm->Bcast(count, rank);
    Array<OneD, char> bytes(count);
    if (comm->GetRank() == rank)
    {
        std::copy(error.begin(), error.end(), bytes.begin());
    }
    comm->Bcast(bytes, rank);
    NEKERROR(ErrorUtil::efatal, std::string(bytes.begin(), bytes.end()));
}
} // namespace

LibUtilities::FieldMetaDataMap HITInflowRestart::Describe(
    const HITInflowConfiguration &configuration,
    const SolverUtils::HITGenerator &generator)
{
    if (!configuration.spectrum || configuration.spectrumDefinition.empty())
    {
        throw std::invalid_argument("HITInflow: incomplete spectrum descriptor");
    }
    const auto &config = generator.GetConfig();
    const auto diagnostics = generator.Diagnostics();
    LibUtilities::FieldMetaDataMap result;
    auto put = [&result](const std::string &key, const std::string &value)
    { result["HITInflow." + key] = value; };
    put("Schema", "1");
    put("Algorithm", SolverUtils::HITGenerator::Algorithm());
    put("Seed", std::to_string(config.seed));
    put("Lengths", Tuple(config.lengths));
    put("Counts", Tuple(config.counts));
    put("Origin", Tuple(config.origin));
    put("Mean", Tuple(SolverUtils::HITVector{config.meanSpeed, 0.0, 0.0}));
    put("TimeOrigin", Number(config.timeOrigin));
    put("Spectrum", configuration.spectrumDefinition);
    put("SpectrumSHA256", configuration.spectrumDigest.empty()
                              ? "none" : configuration.spectrumDigest);
    put("Cutoff", Number(config.kMax));
    put("TaperStart", Number(config.taperStart));
    put("Normalization", config.normalization == SolverUtils::HITNormalization::Preserve
                             ? "Preserve" : "TargetRMS");
    put("TargetRMS", Number(config.targetRMS));
    put("EnergyScale", Number(diagnostics.energyScale));
    put("Environment", Environment());

    std::ostringstream table;
    table.imbue(std::locale::classic());
    const auto &k = configuration.spectrum->GetWavenumbers();
    const auto &energy = configuration.spectrum->GetEnergies();
    for (std::size_t i = 0; i < k.size(); ++i)
    {
        if (i)
        {
            table << ' ';
        }
        table << Number(k[i]) << ' '
              << Number(energy[i] * diagnostics.energyScale);
    }
    put("NormalizedTable", k.empty() ? "none" : table.str());

    std::ostringstream coefficients;
    coefficients.imbue(std::locale::classic());
    for (const auto &mode : generator.GetModes())
    {
        coefficients << Tuple(mode.index) << ' ' << Tuple(mode.wavevector) << ' ';
        for (const auto &value : mode.coefficient)
        {
            coefficients << Number(value.real()) << ' ' << Number(value.imag()) << ' ';
        }
        coefficients << '\n';
    }
    put("ModesSHA256", HITInflowConfig::SHA256(coefficients.str()));
    std::ostringstream identity;
    identity.imbue(std::locale::classic());
    for (const auto &entry : result)
    {
        identity << entry.first.size() << ':' << entry.first
                 << entry.second.size() << ':' << entry.second;
    }
    put("IdentitySHA256", HITInflowConfig::SHA256(identity.str()));
    return result;
}

void HITInflowRestart::Verify(const LibUtilities::FieldMetaDataMap &stored,
                             const LibUtilities::FieldMetaDataMap &expected)
{
    for (const auto &entry : expected)
    {
        const auto found = stored.find(entry.first);
        if (found == stored.end() || found->second != entry.second)
        {
            throw std::invalid_argument(
                "HITInflow: restart metadata missing or mismatched: " + entry.first);
        }
    }
}

void HITInflowRestart::Prepare(
    const LibUtilities::SessionReaderSharedPtr &session,
    const HITInflowConfiguration &configuration,
    const SolverUtils::HITGenerator &generator,
    LibUtilities::FieldMetaDataMap &metadata, NekDouble initialTime)
{
    std::string error;
    LibUtilities::FieldMetaDataMap expected;
    try
    {
        expected = Describe(configuration, generator);
        int fileVariables = 0;
        std::string filename;
        if (session->DefinesFunction("InitialConditions"))
        {
            for (const auto &variable : session->GetVariables())
            {
                if (session->GetFunctionType("InitialConditions", variable) ==
                    LibUtilities::eFunctionTypeFile)
                {
                    const auto input = session->GetFunctionFilename("InitialConditions", variable);
                    if (!filename.empty() && filename != input)
                    {
                        throw std::invalid_argument(
                            "HITInflow: restart u,v,w,p must use one common field file");
                    }
                    filename = input;
                    ++fileVariables;
                }
            }
        }
        if (fileVariables)
        {
            if (fileVariables != 4 || !metadata.count("Time"))
            {
                throw std::invalid_argument(
                    "HITInflow: restart requires all u,v,w,p and stored Time metadata");
            }
            // Native CheckForRestartTime has imported this map already. Do not
            // seed new keys before this verification: native imports merge.
            Verify(metadata, expected);
        }
        const auto &config = generator.GetConfig();
        const long double begin = config.timeOrigin;
        const long double end = begin + config.recordDuration;
        const long double finish = static_cast<long double>(initialTime) +
            session->GetParameter("NumSteps") * static_cast<long double>(config.timeStep);
        const long double tolerance = 64.0L * std::numeric_limits<NekDouble>::epsilon() *
            std::max({1.0L, std::abs(begin), std::abs(end),
                      std::abs(static_cast<long double>(initialTime))});
        if (!std::isfinite(initialTime) || initialTime < begin - tolerance ||
            finish > end + tolerance)
        {
            throw std::invalid_argument(
                "HITInflow: actual job start/end exceeds the declared frozen record");
        }
    }
    catch (const std::exception &exception)
    {
        error = exception.what();
    }
    Agree(session->GetComm(), error);
    for (const auto &entry : expected)
    {
        metadata[entry.first] = entry.second;
    }
    // Informational fields do not define the random realization.
    metadata["HITInflow.BuildGit"] = LibUtilities::GitConsts().GetSha1();
    metadata["HITInflow.RecordDuration"] = Number(generator.GetConfig().recordDuration);
    metadata["HITInflow.Repeat"] = generator.GetConfig().allowRepeat ? "Allow" : "Error";
    metadata["HITInflow.TimeStep"] = Number(generator.GetConfig().timeStep);
    if (session->GetComm()->GetRank() == 0)
    {
        std::cout << "HITInflow restart identity: "
                  << expected.at("HITInflow.IdentitySHA256")
                  << ", start time: " << Number(initialTime) << std::endl;
    }
}
} // namespace Nektar
