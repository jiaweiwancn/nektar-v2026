///////////////////////////////////////////////////////////////////////////////
// File: HITFlowDiagnostics.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#include "HITFieldIO.h"
#include "HITTestFields.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace Nektar;

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        std::cerr << "Usage: HITFlowDiagnostics session.xml solution.fld\n";
        return 2;
    }
    try
    {
        char *arguments[] = {argv[0], argv[1], nullptr};
        auto session = LibUtilities::SessionReader::CreateInstance(2, arguments);
        auto graph = SpatialDomains::MeshGraphIO::Read(session);
        auto fields = HITTestTools::Fields(session, graph);
        auto data = HITTestTools::ImportFields(session, fields, argv[2]);
        const int points = fields[0]->GetNpoints();
        Array<OneD, NekDouble> divergence(points, 0.0), ones(points, 1.0);
        NekDouble velocitySquared = 0.0, gradientSquared = 0.0;
        for (int component = 0; component < 3; ++component)
        {
            Array<OneD, NekDouble> dx(points), dy(points), dz(points);
            fields[component]->PhysDeriv(data.physical[component], dx, dy, dz);
            const auto &diagonal = component == 0 ? dx : component == 1 ? dy : dz;
            for (int q = 0; q < points; ++q)
            {
                divergence[q] += diagonal[q];
            }
            for (const auto &gradient : {dx, dy, dz})
            {
                const auto norm = fields[component]->L2(gradient);
                gradientSquared += norm * norm;
            }
            const auto norm = fields[component]->L2(data.physical[component]);
            velocitySquared += norm * norm;
        }
        const NekDouble volume = fields[0]->Integral(ones);
        const NekDouble velocityNorm = std::sqrt(velocitySquared);
        const NekDouble divNorm = fields[0]->L2(divergence);
        const NekDouble divIntegral = fields[0]->Integral(divergence);
        ASSERTL0(volume > 0.0 && velocityNorm > 0.0,
                 "HIT diagnostic requires a nonzero velocity on a positive volume");
        const NekDouble referenceLength = std::cbrt(volume);
        const NekDouble relativeDivergence = referenceLength * divNorm / velocityNorm;

        // Integrate solution traces, not prescribed boundary-condition arrays.
        // In H1D the z boundary is periodic; integrate each physical xy plane
        // with its uniform Fourier quadrature weight. This tool imports the
        // complete field and runs serially, so no boundary-local collectives
        // are entered by only a subset of MPI ranks.
        const bool homogeneous = fields[0]->GetExpType() == MultiRegions::e3DH1D;
        const int planes = homogeneous ? fields[0]->GetZIDs().size() : 1;
        const NekDouble weight = homogeneous ? fields[0]->GetHomoLen() / planes : 1.0;
        NekDouble inletFlux = 0.0, outletFlux = 0.0;
        int inletRegions = 0, outletRegions = 0;
        for (int plane = 0; plane < planes; ++plane)
        {
            auto geometry = homogeneous ? fields[0]->GetPlane(plane) : fields[0];
            const int planePoints = geometry->GetNpoints();
            for (int region = 0; region < geometry->GetBndConditions().size(); ++region)
            {
                const auto condition = geometry->GetBndConditions()[region];
                if (condition->GetBoundaryConditionType() == SpatialDomains::ePeriodic)
                {
                    continue;
                }
                auto boundary = geometry->GetBndCondExpansions()[region];
                Array<OneD, Array<OneD, NekDouble>> normals;
                geometry->GetBoundaryNormals(region, normals);
                Array<OneD, NekDouble> normalVelocity(boundary->GetNpoints(), 0.0);
                for (int component = 0; component < normals.size(); ++component)
                {
                    Array<OneD, NekDouble> trace;
                    geometry->ExtractPhysToBnd(region,
                        data.physical[component] + plane * planePoints, trace);
                    for (int q = 0; q < trace.size(); ++q)
                    {
                        normalVelocity[q] += trace[q] * normals[component][q];
                    }
                }
                const auto flux = weight * boundary->Integral(normalVelocity);
                if (condition->GetUserDefined() == "HITInflow")
                {
                    inletFlux += flux;
                    ++inletRegions;
                }
                else
                {
                    outletFlux += flux;
                    ++outletRegions;
                }
            }
        }
        ASSERTL0(inletRegions && outletRegions,
                 "HIT diagnostic requires inlet and outlet traces");
        const NekDouble netFlux = inletFlux + outletFlux;
        const NekDouble fluxScale = 0.5 * (std::abs(inletFlux) + std::abs(outletFlux));
        ASSERTL0(fluxScale > 0.0, "HIT diagnostic requires nonzero through-flow");
        const NekDouble relativeFlux = std::abs(netFlux) / fluxScale;
        const NekDouble gauss = std::abs(divIntegral - netFlux) / fluxScale;
        for (auto value : {volume, velocityNorm, divNorm, referenceLength,
                          relativeDivergence, relativeFlux, gauss, gradientSquared})
        {
            ASSERTL0(std::isfinite(value), "Non-finite HIT flow diagnostic");
        }
        std::ofstream output(session->GetSessionName() + ".flow.json");
        output << std::setprecision(17)
               << "{\"volume\":" << volume
               << ",\"velocity_l2\":" << velocityNorm
               << ",\"reference_length\":" << referenceLength
               << ",\"gradient_l2\":" << std::sqrt(gradientSquared)
               << ",\"divergence_l2\":" << divNorm
               << ",\"divergence_integral\":" << divIntegral
               << ",\"relative_divergence\":" << relativeDivergence
               << ",\"inlet_flux\":" << inletFlux
               << ",\"outlet_flux\":" << outletFlux
               << ",\"relative_flux_imbalance\":" << relativeFlux
               << ",\"gauss_relative_residual\":" << gauss << "}\n";
        ASSERTL0(output.good(), "Cannot write HIT flow diagnostic");
        std::cout << std::setprecision(17)
                  << "HIT_FLOW div_l2 " << divNorm
                  << " div_relative " << relativeDivergence
                  << " flux_relative " << relativeFlux
                  << " gauss_relative " << gauss << '\n';
        session->Finalise();
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << std::endl;
        return 1;
    }
}
