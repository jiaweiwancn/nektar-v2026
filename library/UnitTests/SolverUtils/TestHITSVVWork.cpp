///////////////////////////////////////////////////////////////////////////////
// File: TestHITSVVWork.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#include "HITSVVWork.h"
#include "HITTestFields.h"
#include <LibUtilities/BasicUtils/Filesystem.hpp>
#include <LibUtilities/Communication/CommSerial.h>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>
#include <boost/test/unit_test.hpp>
#include <tinyxml.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace Nektar::HITSVVWorkTests
{
namespace
{
constexpr NekDouble amplitude = 0.05;
constexpr NekDouble wave = 3.0;
using Gradients = std::array<std::array<Array<OneD, NekDouble>, 3>, 3>;

struct TemporaryDirectory
{
    fs::path path = fs::temp_directory_path() /
                    LibUtilities::UniquePath("hit-svv-work");
    TemporaryDirectory()
    {
        if (!fs::create_directory(path))
        {
            throw std::runtime_error("Cannot create SVV work test directory");
        }
    }
    ~TemporaryDirectory()
    {
        try
        {
            fs::remove_all(path);
        }
        catch (...)
        {
        }
    }
};

struct Fixture
{
    TemporaryDirectory directory;
    LibUtilities::SessionReaderSharedPtr session;
    Array<OneD, MultiRegions::ExpListSharedPtr> fields;
    NekDouble volume = 0.0;

    Fixture()
    {
        const auto source = fs::path(HIT_TEST_FIXTURE_DIR) /
                            "m6-svv-fourier-pair.xml";
        TiXmlDocument document;
        BOOST_REQUIRE(document.LoadFile(source.string().c_str()));
        auto conditions = document.RootElement()->FirstChildElement("CONDITIONS");
        BOOST_REQUIRE(conditions);
        auto info = conditions->FirstChildElement("SOLVERINFO");
        BOOST_REQUIRE(info);
        TiXmlElement dg("I");
        dg.SetAttribute("PROPERTY", "SpectralVanishingViscositySpectralHP");
        dg.SetAttribute("VALUE", "DGKernel");
        info->InsertEndChild(dg);
        auto parameters = conditions->FirstChildElement("PARAMETERS");
        BOOST_REQUIRE(parameters);
        TiXmlElement coefficient("P");
        coefficient.LinkEndChild(new TiXmlText("SVVDiffCoeff = 1"));
        parameters->InsertEndChild(coefficient);
        const auto input = directory.path / "probe.xml";
        BOOST_REQUIRE(document.SaveFile(input.string().c_str()));

        // Collections may write tuning files beside a session. All inputs and
        // generated tuning files belong to the temporary directory.
        char program[] = "HITSVVWorkUnitTests";
        char *arguments[] = {program, nullptr};
        std::vector<std::string> inputs{input.string()};
        auto comm = LibUtilities::CommSerial::create(1, arguments);
        session = LibUtilities::SessionReader::CreateInstance(
            1, arguments, inputs, comm);
        auto graph = SpatialDomains::MeshGraphIO::Read(session);
        fields = HITTestTools::Fields(session, graph);
        BOOST_REQUIRE_EQUAL(fields[0]->GetZIDs().size(), 8U);
        Array<OneD, NekDouble> ones(fields[0]->GetNpoints(), 1.0);
        volume = fields[0]->Integral(ones);
        BOOST_REQUIRE(std::isfinite(volume) && volume > 0.0);
    }

    Gradients Sample(bool polynomial, bool paired, NekDouble phase = 0.0)
    {
        const int count = fields[0]->GetNpoints();
        Array<OneD, NekDouble> x(count), y(count), z(count);
        fields[0]->GetCoords(x, y, z);
        for (int point = 0; point < count; ++point)
        {
            fields[0]->UpdatePhys()[point] =
                (polynomial ? x[point] * x[point] * x[point] : 0.0) +
                amplitude * std::cos(wave * z[point] + phase);
            fields[1]->UpdatePhys()[point] =
                (polynomial ? y[point] * y[point] * y[point] : 0.0) +
                (paired ? amplitude * std::sin(wave * z[point] + phase) : 0.0);
            fields[2]->UpdatePhys()[point] = 0.0;
            fields[3]->UpdatePhys()[point] = 0.0;
        }
        Gradients gradients;
        for (int component = 0; component < 3; ++component)
        {
            fields[component]->SetWaveSpace(false);
            fields[component]->SetPhysState(true);
            for (int direction = 0; direction < 3; ++direction)
            {
                gradients[component][direction] = Array<OneD, NekDouble>(count);
            }
            fields[component]->PhysDeriv(fields[component]->GetPhys(),
                gradients[component][0], gradients[component][1],
                gradients[component][2]);
        }
        return gradients;
    }

    NekDouble MatrixWork() const
    {
        const int planes = fields[0]->GetZIDs().size();
        const int planePoints = fields[0]->GetPlane(0)->GetNpoints();
        const NekDouble weight = fields[0]->GetHomoLen() / planes / volume;
        long double result = 0.0L;
        for (int plane = 0; plane < planes; ++plane)
        {
            auto geometry = fields[0]->GetPlane(plane);
            for (int element = 0; element < geometry->GetExpSize(); ++element)
            {
                auto expansion = geometry->GetExp(element);
                NekDouble h = 0.0;
                for (int edge = 0; edge < expansion->GetGeom()->GetNumEdges(); ++edge)
                {
                    const auto segment = expansion->GetGeom()->GetEdge(edge);
                    h = std::max(h, segment->GetVertex(0)->dist(
                        *segment->GetVertex(1)));
                }
                const int degree = std::max(expansion->GetBasisNumModes(0),
                    expansion->GetBasisNumModes(1)) - 1;
                BOOST_REQUIRE_GT(degree, 0);
                StdRegions::ConstFactorMap factors;
                factors[StdRegions::eFactorSVVCutoffRatio] = 0.75;
                factors[StdRegions::eFactorSVVDiffCoeff] = 1.0;
                // Default native DG velocity magnitude is one. The variable
                // diffusivity contributes h/p; the physical coefficient is 1.
                factors[StdRegions::eFactorSVVDGKerDiffCoeff] = h / degree;
                const auto augmented = expansion->GetLocMatrix(
                    StdRegions::eLaplacian, factors);
                const auto baseline = expansion->GetLocMatrix(
                    StdRegions::eLaplacian);
                const int coefficients = expansion->GetNcoeffs();
                const int points = expansion->GetTotPoints();
                const int offset = plane * planePoints +
                                   geometry->GetPhys_Offset(element);
                BOOST_REQUIRE_EQUAL(augmented->GetRows(), coefficients);
                BOOST_REQUIRE_EQUAL(baseline->GetRows(), coefficients);
                for (int component = 0; component < 3; ++component)
                {
                    Array<OneD, NekDouble> physical(points), modal(coefficients);
                    for (int q = 0; q < points; ++q)
                    {
                        physical[q] = fields[component]->GetPhys()[offset + q];
                    }
                    // Local polynomial operator probe: no global periodic or
                    // Dirichlet projection changes these sampled polynomials.
                    expansion->FwdTrans(physical, modal);
                    long double quadratic = 0.0L;
                    for (int row = 0; row < coefficients; ++row)
                    {
                        for (int column = 0; column < coefficients; ++column)
                        {
                            quadratic += static_cast<long double>(modal[row]) *
                                ((*augmented)(row, column) -
                                 (*baseline)(row, column)) * modal[column];
                        }
                    }
                    result += weight * quadratic;
                }
            }
        }
        return static_cast<NekDouble>(result);
    }
};
} // namespace

BOOST_AUTO_TEST_SUITE(HITSVVWork)

BOOST_AUTO_TEST_CASE(DGWorkMatchesNativeStiffnessQuadraticForm)
{
    Fixture fixture;
    const auto gradients = fixture.Sample(true, true);
    const auto work = HITTestTools::SVVWork(fixture.session, fixture.fields[0],
                                           gradients, fixture.volume);
    const auto matrixWork = fixture.MatrixWork();
    const NekDouble expectedFourier = amplitude * amplitude * wave * wave *
                                      0.1 * std::exp(-1.0);
    BOOST_REQUIRE_GT(matrixWork, 0.0);
    BOOST_REQUIRE_GT(work[0], 0.0);
    BOOST_REQUIRE_GT(work[1], 0.0);
    BOOST_CHECK_CLOSE_FRACTION(work[0], matrixWork, 1e-9);
    BOOST_CHECK_CLOSE_FRACTION(work[1], expectedFourier, 1e-9);
    BOOST_TEST_MESSAGE("SVV xy bulk work: helper=" << work[0]
        << ", native stiffness difference=" << matrixWork
        << "; Fourier work=" << work[1] << ", analytic=" << expectedFourier);
}

BOOST_AUTO_TEST_CASE(FourierWorkAndParsevalAreIndependentOfPhase)
{
    Fixture fixture;
    const NekDouble variance = 0.5 * amplitude * amplitude * wave * wave;
    const NekDouble expectedWork = 0.1 * std::exp(-1.0) * variance;
    for (const NekDouble phase : {0.0, std::acos(-1.0) / 2.0, 0.37})
    {
        BOOST_TEST_CONTEXT("phase=" << phase)
        {
            const auto gradients = fixture.Sample(false, false, phase);
            const auto work = HITTestTools::SVVWork(fixture.session,
                fixture.fields[0], gradients, fixture.volume);
            const auto derivativeNorm = fixture.fields[0]->L2(gradients[0][2]);
            BOOST_CHECK_CLOSE_FRACTION(derivativeNorm * derivativeNorm /
                                       fixture.volume, variance, 1e-10);
            BOOST_CHECK_CLOSE_FRACTION(work[1], expectedWork, 1e-10);
            BOOST_CHECK_SMALL(work[0], 1e-20);
        }
    }
}

BOOST_AUTO_TEST_SUITE_END()
} // namespace Nektar::HITSVVWorkTests
