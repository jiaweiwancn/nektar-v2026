///////////////////////////////////////////////////////////////////////////////
// File: TestHITBoundary.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/BoundaryConditions/IncBoundaryConditions.h>
#include <LibUtilities/Communication/CommSerial.h>
#include <MultiRegions/ContField.h>
#include <SpatialDomains/MeshGraphIO.h>
#include <boost/test/unit_test.hpp>
#include <tinyxml.h>
#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>

namespace Nektar::HITBoundaryTests
{
namespace
{
LibUtilities::SessionReaderSharedPtr Session(const std::string &fixture)
{
    char name[] = "HITBoundaryUnitTests";
    char *arguments[] = {name, nullptr};
    std::vector<std::string> files{std::string(HIT_TEST_FIXTURE_DIR) + "/" +
                                   fixture};
    auto comm = LibUtilities::CommSerial::create(1, arguments);
    auto session = LibUtilities::SessionReader::CreateInstance(1, arguments,
                                                               files, comm);
    session->InitSession();
    return session;
}

TiXmlElement *Region(const LibUtilities::SessionReaderSharedPtr &session)
{
    return session->GetDocument().RootElement()->FirstChildElement("HITINFLOW")
        ->FirstChildElement("REGION");
}

struct BoundaryFixture
{
    LibUtilities::SessionReaderSharedPtr session;
    Array<OneD, MultiRegions::ExpListSharedPtr> fields;
    IncBoundaryConditions manager;
    std::shared_ptr<const HITInflowConfiguration> configuration;
    std::shared_ptr<SolverUtils::HITGenerator> reference;
    int inlet = -1;

    explicit BoundaryFixture(
        const std::string &file,
        const std::function<void(TiXmlElement *)> &change = {})
        : session(Session(file)), fields(4)
    {
        // Resolve the smooth periodic modes accurately enough to measure flux
        // independently of low-degree spatial truncation.
        auto graph = SpatialDomains::MeshGraphIO::Read(session);
        graph->ReadExpansionInfo();
        graph->SetExpansionInfoToNumModes(12);
        for (int component = 0; component < 4; ++component)
        {
            fields[component] =
                MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(
                    session, graph, session->GetVariable(component), true,
                    true);
        }
        if (change)
        {
            change(Region(session));
        }
        configuration = HITInflowConfig::Load(session);
        reference = std::make_shared<SolverUtils::HITGenerator>(
            configuration->generator, *configuration->spectrum);
        manager.Initialize(session, fields);
        for (int boundary = 0;
             boundary < fields[0]->GetBndConditions().size(); ++boundary)
        {
            if (fields[0]->GetBndConditions()[boundary]->GetUserDefined() ==
                "HITInflow")
            {
                inlet = boundary;
            }
        }
        BOOST_REQUIRE(inlet >= 0);
        BOOST_REQUIRE_EQUAL(configuration->regionId, 7);
        BOOST_REQUIRE(inlet != configuration->regionId);
    }

    std::vector<NekDouble> Coefficients(int component) const
    {
        const auto &coefficients = fields[component]
                                       ->GetBndCondExpansions()[inlet]
                                       ->GetCoeffs();
        return {coefficients.begin(), coefficients.end()};
    }

    void CheckSamples(NekDouble time)
    {
        auto boundary = fields[0]->GetBndCondExpansions()[inlet];
        const int points = boundary->GetTotPoints();
        Array<OneD, NekDouble> x(points), y(points), z(points);
        boundary->GetCoords(x, y, z);
        for (int point = 0; point < points; ++point)
        {
            const auto expected = reference->Evaluate(
                {x[point], y[point], z[point]}, time);
            for (int component = 0; component < 3; ++component)
            {
                BOOST_CHECK_EQUAL(
                    fields[component]->GetBndCondExpansions()[inlet]
                        ->GetPhys()[point],
                    expected.velocity[component]);
            }
        }
    }
};
} // namespace

BOOST_AUTO_TEST_SUITE(HITBoundary)

BOOST_AUTO_TEST_CASE(SHA256KnownAnswers)
{
    BOOST_CHECK_EQUAL(HITInflowConfig::SHA256(""),
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
    BOOST_CHECK_EQUAL(HITInflowConfig::SHA256("abc"),
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    BOOST_CHECK_EQUAL(HITInflowConfig::SHA256(std::string(1000000, 'a')),
        "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0");
}

BOOST_AUTO_TEST_CASE(StrictConfigurationFailures)
{
    const std::vector<std::function<void(TiXmlElement *)>> changes{
        [](auto root) { root->SetAttribute("UNKNOWN", "1"); },
        [](auto root) { root->FirstChildElement("SYNTHESIS")
                           ->SetAttribute("SEED", "18446744073709551616"); },
        [](auto root) { root->FirstChildElement("SYNTHESIS")
                           ->SetAttribute("NX", "5"); },
        [](auto root) { root->FirstChildElement("ORIGIN")
                           ->SetAttribute("X", "t"); },
        [](auto root) { root->FirstChildElement("SPECTRUM")
                           ->SetAttribute("AMPLITUDE", "1"); },
        [](auto root) { root->FirstChildElement("CUTOFF")
                           ->SetAttribute("KMAX", "100"); },
        [](auto root) { root->InsertEndChild(*root->FirstChildElement("ORIGIN")); }
    };
    for (const auto &change : changes)
    {
        auto session = Session("m3-mean-only.xml");
        change(Region(session));
        BOOST_CHECK_THROW(HITInflowConfig::Load(session), std::exception);
    }
    auto session = Session("m3-mean-only.xml");
    auto conditions = session->GetDocument().RootElement()
                          ->FirstChildElement("CONDITIONS");
    conditions->FirstChildElement("BOUNDARYCONDITIONS")
        ->FirstChildElement("REGION")->FirstChildElement("D")
        ->SetAttribute("VALUE", "1+t");
    BOOST_CHECK_THROW(HITInflowConfig::Load(session), std::exception);
}

BOOST_AUTO_TEST_CASE(MeanOnlyProjectedBoundary)
{
    BoundaryFixture fixture("m3-mean-only.xml");
    BOOST_CHECK(fixture.reference->GetModes().empty());
    fixture.manager.UpdatePrescribedVelocity(0.0);
    fixture.CheckSamples(0.0);
    for (int component = 0; component < 3; ++component)
    {
        auto boundary = fixture.fields[component]
                            ->GetBndCondExpansions()[fixture.inlet];
        Array<OneD, NekDouble> projected(boundary->GetTotPoints());
        boundary->BwdTrans(boundary->GetCoeffs(), projected);
        for (NekDouble value : projected)
        {
            BOOST_CHECK_SMALL(value - (component == 0 ? 1.0 : 0.0), 2e-12);
        }
    }
    fixture.manager.UpdatePrescribedVelocity(0.02);
    fixture.CheckSamples(0.02);
}

BOOST_AUTO_TEST_CASE(RejectInconsistentGeometry)
{
    BOOST_CHECK_THROW(
        BoundaryFixture("m3-mean-only.xml", [](auto region) {
            region->FirstChildElement("ORIGIN")->SetAttribute("X", "0.1");
        }), std::exception);
    BOOST_CHECK_THROW(
        BoundaryFixture("m3-mean-only.xml", [](auto region) {
            region->FirstChildElement("SYNTHESIS")->SetAttribute("LY", "6");
        }), std::exception);
}

BOOST_AUTO_TEST_CASE(NoHITConfigurationIsOptional)
{
    auto session = Session("m3-mean-only.xml");
    auto root = session->GetDocument().RootElement();
    root->RemoveChild(root->FirstChildElement("HITINFLOW"));
    auto inlet = root->FirstChildElement("CONDITIONS")
                     ->FirstChildElement("BOUNDARYCONDITIONS")
                     ->FirstChildElement("REGION");
    for (auto velocity = inlet->FirstChildElement("D"); velocity;
         velocity = velocity->NextSiblingElement("D"))
    {
        velocity->RemoveAttribute("USERDEFINEDTYPE");
    }
    BOOST_CHECK(!HITInflowConfig::Load(session));
}

BOOST_AUTO_TEST_CASE(ThreeComponentsTimeAndRepeatability)
{
    BoundaryFixture fixture("m3-stochastic.xml");
    const auto pressureBefore = fixture.Coefficients(3);
    fixture.manager.UpdatePrescribedVelocity(0.0);
    fixture.CheckSamples(0.0);
    const std::vector<std::vector<NekDouble>> before{
        fixture.Coefficients(0), fixture.Coefficients(1),
        fixture.Coefficients(2)};
    fixture.manager.UpdatePrescribedVelocity(0.01);
    fixture.CheckSamples(0.01);
    for (int component = 0; component < 3; ++component)
    {
        BOOST_CHECK(fixture.Coefficients(component) != before[component]);
    }
    fixture.manager.UpdatePrescribedVelocity(0.02);
    fixture.CheckSamples(0.02);
    fixture.manager.UpdatePrescribedVelocity(0.0);
    for (int component = 0; component < 3; ++component)
    {
        BOOST_CHECK(fixture.Coefficients(component) == before[component]);
    }
    BOOST_CHECK(fixture.Coefficients(3) == pressureBefore);

    auto boundary = fixture.fields[0]->GetBndCondExpansions()[fixture.inlet];
    Array<OneD, NekDouble> fluctuation(boundary->GetTotPoints());
    boundary->BwdTrans(boundary->GetCoeffs(), fluctuation);
    for (auto &value : fluctuation)
    {
        value -= 1.0;
    }
    const auto flux = boundary->Integral(fluctuation);
    BOOST_TEST_MESSAGE("Projected inlet fluctuating flux = " << flux);
    BOOST_CHECK_SMALL(flux, 2e-8);
}

BOOST_AUTO_TEST_CASE(PrescribedNormalAcceleration)
{
    BoundaryFixture fixture("m3-stochastic.xml");
    constexpr NekDouble time = 0.01, step = 1e-6;
    auto boundary = fixture.fields[0]->GetBndCondExpansions()[fixture.inlet];
    const int points = boundary->GetTotPoints();
    Array<OneD, NekDouble> x(points), y(points), z(points);
    boundary->GetCoords(x, y, z);
    fixture.manager.UpdatePrescribedVelocity(time - step);
    std::vector<NekDouble> previous(boundary->GetPhys().begin(),
                                    boundary->GetPhys().end());
    fixture.manager.UpdatePrescribedVelocity(time + step);
    NekDouble maximum = 0.0;
    for (int point = 0; point < points; ++point)
    {
        const NekDouble derivative = fixture.reference
                                        ->Evaluate({x[point], y[point], z[point]},
                                                   time)
                                        .timeDerivative[0];
        maximum = std::max(maximum, std::abs(derivative));
        const NekDouble difference =
            (boundary->GetPhys()[point] - previous[point]) / (2.0 * step);
        BOOST_CHECK_SMALL(difference - derivative, 2e-9);
    }
    BOOST_CHECK_GT(maximum, 1e-3);
}

BOOST_AUTO_TEST_SUITE_END()
} // namespace Nektar::HITBoundaryTests
