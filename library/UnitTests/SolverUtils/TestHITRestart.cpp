///////////////////////////////////////////////////////////////////////////////
// File: TestHITRestart.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/BoundaryConditions/HITInflowRestart.h>
#include <boost/test/unit_test.hpp>
#include <tinyxml.h>

#include <algorithm>
#include <functional>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Nektar::HITRestartTests
{
namespace
{
using SolverUtils::HITGenerator;
using SolverUtils::HITHighWavenumber;
using SolverUtils::HITInterpolation;
using SolverUtils::HITLowWavenumber;
using SolverUtils::HITNormalization;
using SolverUtils::HITSpectrum;
using Metadata = LibUtilities::FieldMetaDataMap;

const std::string prefix = "HITInflow.";
const std::string tableHeader =
    "# NEKTAR_HIT_SPECTRUM_V1 KIND=Shell3D K_UNITS=rad/m E_UNITS=m3/s2\n";
const std::string tableRows = "0.5 0.02\n1 0.1\n1.5 0.2\n2 0\n";

HITInflowConfiguration VonKarmanConfiguration()
{
    HITInflowConfiguration configuration;
    configuration.regionId = 7;
    configuration.generator.counts = {4, 4, 4};
    configuration.generator.seed = 20261028;
    configuration.generator.kMax = 1.55;
    configuration.generator.taperStart = 1.5;
    configuration.generator.normalization = HITNormalization::TargetRMS;
    configuration.generator.targetRMS = 0.05;
    configuration.spectrum = std::make_shared<HITSpectrum>(
        HITSpectrum::VonKarman(0.8, 1.0));
    configuration.spectrumDefinition = "VonKarman LENGTH=0.8 AMPLITUDE=1";
    configuration.description = "Test job one";
    return configuration;
}

HITInflowConfiguration TabulatedConfiguration(const std::string &rows = tableRows)
{
    auto configuration = VonKarmanConfiguration();
    const std::string bytes = tableHeader + rows;
    std::istringstream input(bytes);
    configuration.spectrum = std::make_shared<HITSpectrum>(
        HITSpectrum::ReadTable(input, HITInterpolation::Linear,
                               HITLowWavenumber::Zero, HITHighWavenumber::Zero));
    configuration.spectrumDigest = HITInflowConfig::SHA256(bytes);
    configuration.spectrumDefinition =
        "Tabulated INTERPOLATION=Linear LOW=Zero HIGH=Zero";
    return configuration;
}

Metadata Describe(const HITInflowConfiguration &configuration)
{
    const HITGenerator generator(configuration.generator,
                                  *configuration.spectrum);
    return HITInflowRestart::Describe(configuration, generator);
}

bool IsSHA256(const std::string &text)
{
    return text.size() == 64 &&
           std::all_of(text.begin(), text.end(), [](char character) {
               return (character >= '0' && character <= '9') ||
                      (character >= 'a' && character <= 'f');
           });
}

void CheckDifferentIdentity(const Metadata &original, const Metadata &changed)
{
    BOOST_CHECK(original != changed);
    BOOST_CHECK(original.at(prefix + "IdentitySHA256") !=
                changed.at(prefix + "IdentitySHA256"));
    BOOST_CHECK_THROW(HITInflowRestart::Verify(original, changed),
                      std::invalid_argument);
    BOOST_CHECK_THROW(HITInflowRestart::Verify(changed, original),
                      std::invalid_argument);
}
} // namespace

BOOST_AUTO_TEST_SUITE(HITRestart)

BOOST_AUTO_TEST_CASE(CompleteDeterministicDescription)
{
    const auto configuration = VonKarmanConfiguration();
    const auto first = Describe(configuration);
    const auto second = Describe(configuration);
    BOOST_CHECK(first == second);
    BOOST_CHECK_NO_THROW(HITInflowRestart::Verify(first, second));

    for (const std::string key : {
             "Schema", "Algorithm", "Seed", "Lengths", "Counts", "Origin",
             "Mean", "TimeOrigin", "Spectrum", "SpectrumSHA256",
             "NormalizedTable", "Cutoff", "TaperStart", "Normalization",
             "TargetRMS", "EnergyScale", "Environment", "ModesSHA256",
             "IdentitySHA256"})
    {
        BOOST_TEST_CONTEXT(key)
        {
            BOOST_CHECK_EQUAL(first.count(prefix + key), 1U);
        }
    }
    BOOST_CHECK_EQUAL(first.at(prefix + "Algorithm"), HITGenerator::Algorithm());
    BOOST_CHECK_EQUAL(first.at(prefix + "Seed"), "20261028");
    BOOST_CHECK_EQUAL(first.at(prefix + "NormalizedTable"), "none");
    BOOST_CHECK(IsSHA256(first.at(prefix + "ModesSHA256")));
    BOOST_CHECK(IsSHA256(first.at(prefix + "IdentitySHA256")));
    BOOST_CHECK(!first.at(prefix + "Environment").empty());
}

BOOST_AUTO_TEST_CASE(EveryRequiredKeyIsVerified)
{
    const auto expected = Describe(VonKarmanConfiguration());
    BOOST_REQUIRE(!expected.empty());
    for (const auto &entry : expected)
    {
        BOOST_TEST_CONTEXT(entry.first)
        {
            auto stored = expected;
            stored.erase(entry.first);
            BOOST_CHECK_THROW(HITInflowRestart::Verify(stored, expected),
                              std::invalid_argument);
            stored = expected;
            stored[entry.first] = "tampered:" + entry.second;
            BOOST_CHECK_THROW(HITInflowRestart::Verify(stored, expected),
                              std::invalid_argument);
        }
    }
    BOOST_CHECK_THROW(HITInflowRestart::Verify(Metadata{}, expected),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(RealizationParametersChangeIdentity)
{
    const auto configuration = VonKarmanConfiguration();
    const auto original = Describe(configuration);
    using Change = std::pair<std::string,
                             std::function<void(HITInflowConfiguration &)>>;
    const std::vector<Change> changes{
        {"seed", [](auto &value) { ++value.generator.seed; }},
        {"mean speed", [](auto &value) { value.generator.meanSpeed = 1.25; }},
        {"origin x", [](auto &value) { value.generator.origin[0] = 0.1; }},
        {"origin y", [](auto &value) { value.generator.origin[1] = 0.2; }},
        {"origin z", [](auto &value) { value.generator.origin[2] = 0.3; }},
        {"period x", [](auto &value) { value.generator.lengths[0] += 0.2; }},
        {"period y", [](auto &value) { value.generator.lengths[1] += 0.2; }},
        {"period z", [](auto &value) { value.generator.lengths[2] += 0.2; }},
        {"count x", [](auto &value) { value.generator.counts[0] += 2; }},
        {"count y", [](auto &value) { value.generator.counts[1] += 2; }},
        {"count z", [](auto &value) { value.generator.counts[2] += 2; }},
        {"cutoff", [](auto &value) { value.generator.kMax = 1.9; }},
        {"taper", [](auto &value) { value.generator.taperStart = 1.2; }},
        {"phase time", [](auto &value) { value.generator.timeOrigin = 0.025; }},
        {"normalization", [](auto &value) {
             value.generator.normalization = HITNormalization::Preserve;
         }},
        {"target RMS", [](auto &value) { value.generator.targetRMS = 0.075; }},
        {"spectrum length", [](auto &value) {
             value.spectrum = std::make_shared<HITSpectrum>(
                 HITSpectrum::VonKarman(0.9, 1.0));
             value.spectrumDefinition = "VonKarman LENGTH=0.9 AMPLITUDE=1";
         }}
    };
    for (const auto &change : changes)
    {
        BOOST_TEST_CONTEXT(change.first)
        {
            auto changed = configuration;
            change.second(changed);
            CheckDifferentIdentity(original, Describe(changed));
        }
    }
}

BOOST_AUTO_TEST_CASE(JobSettingsDoNotChangeRealization)
{
    const auto configuration = VonKarmanConfiguration();
    const auto original = Describe(configuration);
    auto changed = configuration;
    changed.generator.recordDuration = 0.75;
    changed.generator.allowRepeat = true;
    changed.generator.timeStep = 0.002;
    changed.generator.maxModeBytes = 1024 * 1024;
    changed.regionId = 123;
    changed.description = "Different record, output schedule and partition";
    const auto current = Describe(changed);
    BOOST_CHECK(original == current);
    BOOST_CHECK_NO_THROW(HITInflowRestart::Verify(original, current));

    auto stored = original;
    stored["Time"] = "0.125";
    stored["ChkFileNum"] = "4";
    stored["SessionName0"] = "previous-session.xml";
    stored["HITInflow.BuildGit"] = "old build provenance";
    stored["HITInflow.RecordDuration"] = "0.125";
    stored["HITInflow.Repeat"] = "Allow";
    stored["ArbitraryProvenance"] = "kept independently of required identity";
    BOOST_CHECK_NO_THROW(HITInflowRestart::Verify(stored, current));
}

BOOST_AUTO_TEST_CASE(NormalizedTableRecordsNumericSpectrum)
{
    const auto configuration = TabulatedConfiguration();
    const HITGenerator generator(configuration.generator,
                                  *configuration.spectrum);
    const auto metadata = HITInflowRestart::Describe(configuration, generator);
    const auto scale = generator.Diagnostics().energyScale;
    BOOST_CHECK_EQUAL(std::stod(metadata.at(prefix + "EnergyScale")), scale);
    BOOST_CHECK_EQUAL(metadata.at(prefix + "SpectrumSHA256"),
                      configuration.spectrumDigest);

    std::istringstream normalized(metadata.at(prefix + "NormalizedTable"));
    const auto &wavenumbers = configuration.spectrum->GetWavenumbers();
    const auto &energies = configuration.spectrum->GetEnergies();
    for (std::size_t row = 0; row < wavenumbers.size(); ++row)
    {
        NekDouble wavenumber, energy;
        BOOST_REQUIRE(static_cast<bool>(normalized >> wavenumber >> energy));
        BOOST_CHECK_EQUAL(wavenumber, wavenumbers[row]);
        BOOST_CHECK_EQUAL(energy, energies[row] * scale);
    }
    std::string extra;
    BOOST_CHECK(!(normalized >> extra));
    BOOST_CHECK_NO_THROW(HITInflowRestart::Verify(metadata, Describe(configuration)));
}

BOOST_AUTO_TEST_CASE(TableValuesAndExactInputBytesAreVerified)
{
    const auto configuration = TabulatedConfiguration();
    const auto original = Describe(configuration);
    const auto changedValues =
        TabulatedConfiguration("0.5 0.02\n1 0.13\n1.5 0.2\n2 0\n");
    const auto changedMetadata = Describe(changedValues);
    CheckDifferentIdentity(original, changedMetadata);
    BOOST_CHECK(original.at(prefix + "NormalizedTable") !=
                changedMetadata.at(prefix + "NormalizedTable"));
    BOOST_CHECK(original.at(prefix + "ModesSHA256") !=
                changedMetadata.at(prefix + "ModesSHA256"));

    auto changedBytes = configuration;
    changedBytes.spectrumDigest = HITInflowConfig::SHA256(
        tableHeader + tableRows + "# Same numbers, changed source bytes.\n");
    const auto byteMetadata = Describe(changedBytes);
    CheckDifferentIdentity(original, byteMetadata);
    BOOST_CHECK_EQUAL(original.at(prefix + "NormalizedTable"),
                      byteMetadata.at(prefix + "NormalizedTable"));
    BOOST_CHECK_EQUAL(original.at(prefix + "ModesSHA256"),
                      byteMetadata.at(prefix + "ModesSHA256"));
}

BOOST_AUTO_TEST_CASE(InvalidSpectrumDescriptorsFail)
{
    auto configuration = VonKarmanConfiguration();
    const HITGenerator generator(configuration.generator,
                                  *configuration.spectrum);
    configuration.spectrumDefinition.clear();
    BOOST_CHECK_THROW(HITInflowRestart::Describe(configuration, generator),
                      std::invalid_argument);
    configuration = VonKarmanConfiguration();
    configuration.spectrum.reset();
    BOOST_CHECK_THROW(HITInflowRestart::Describe(configuration, generator),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(MetadataSurvivesXMLWhitespaceProcessing)
{
    const auto expected = Describe(TabulatedConfiguration());
    TiXmlDocument document;
    auto root = new TiXmlElement("NEKTAR");
    document.LinkEndChild(root);
    auto metadata = new TiXmlElement("Metadata");
    root->LinkEndChild(metadata);
    for (const auto &entry : expected)
    {
        auto element = new TiXmlElement(entry.first);
        element->LinkEndChild(new TiXmlText(entry.second));
        metadata->LinkEndChild(element);
    }
    TiXmlPrinter printer;
    document.Accept(&printer);
    TiXmlDocument restored;
    restored.Parse(printer.CStr());
    BOOST_REQUIRE(!restored.Error());
    Metadata values;
    for (auto node = restored.RootElement()->FirstChildElement("Metadata")
                         ->FirstChildElement(); node;
         node = node->NextSiblingElement())
    {
        BOOST_REQUIRE(node->GetText());
        values[node->ValueStr()] = node->GetText();
    }
    BOOST_CHECK_NO_THROW(HITInflowRestart::Verify(values, expected));
}

BOOST_AUTO_TEST_SUITE_END()
} // namespace Nektar::HITRestartTests
