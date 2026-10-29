///////////////////////////////////////////////////////////////////////////////
// File: HITFlowDiagnostics.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#include "HITFieldIO.h"
#include "HITTestFields.h"
#include "HITSVVWork.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <utility>
#include <vector>

using namespace Nektar;

namespace
{
void WriteSlabs(
    const LibUtilities::SessionReaderSharedPtr &session,
    const MultiRegions::ExpListSharedPtr &field,
    const Array<OneD, Array<OneD, NekDouble>> &physical,
    const Array<OneD, NekDouble> &strainSquared,
    const Array<OneD, NekDouble> &divergence,
    const std::array<Array<OneD, NekDouble>, 2> &svvDensity,
    NekDouble viscosity, NekDouble volume)
{
    ASSERTL0(session->GetComm()->GetSize() == 1 &&
                 field->GetExpType() == MultiRegions::e3DH1D,
             "Slab profiles require a complete serial H1D field");
    NekDouble mean;
    session->LoadParameter("Ubar", mean);
    ASSERTL0(std::isfinite(mean), "Slab profiles require finite prescribed Ubar");
    const auto geometry = field->GetPlane(0);
    const int planes = field->GetZIDs().size();
    const int planePoints = geometry->GetNpoints();
    const NekDouble weight = field->GetHomoLen() / planes;
    using Interval = std::pair<NekDouble, NekDouble>;
    std::vector<Interval> intervals(geometry->GetExpSize());
    std::map<Interval, std::array<NekDouble, 9>> slabs;
    const NekDouble large = std::numeric_limits<NekDouble>::max();
    NekDouble yMin = large, yMax = -large;
    for (int element = 0; element < geometry->GetExpSize(); ++element)
    {
        const auto expansion = geometry->GetExp(element);
        const auto shape = expansion->GetGeom();
        ASSERTL0(expansion->DetShapeType() == LibUtilities::eQuadrilateral &&
                     expansion->GetGeomFactors()->GetGtype() == SpatialDomains::eRegular,
                 "Slab profiles require affine quadrilaterals");
        std::array<NekDouble, 2> lower{large, large}, upper{-large, -large};
        for (int vertex = 0; vertex < 4; ++vertex)
        {
            std::array<NekDouble, 3> coordinates;
            shape->GetVertex(vertex)->GetCoords(coordinates[0], coordinates[1], coordinates[2]);
            for (int d = 0; d < 2; ++d)
            {
                lower[d] = std::min(lower[d], coordinates[d]);
                upper[d] = std::max(upper[d], coordinates[d]);
            }
        }
        ASSERTL0(upper[0] > lower[0] && upper[1] > lower[1],
                 "Degenerate slab element");
        // An affine quadrilateral fills its coordinate bounding rectangle
        // exactly only when its edges align with the x and y axes.
        Array<OneD, NekDouble> localOnes(expansion->GetTotPoints(), 1.0);
        const auto area = expansion->Integral(localOnes);
        const auto rectangle = (upper[0] - lower[0]) * (upper[1] - lower[1]);
        ASSERTL0(area > 0 && std::abs(area - rectangle) <= 1e-10 * rectangle,
                 "Slab profiles require axis-aligned quadrilaterals");
        intervals[element] = {lower[0], upper[0]};
        slabs[intervals[element]] = {};
        yMin = std::min(yMin, lower[1]);
        yMax = std::max(yMax, upper[1]);
    }
    // Integrate complete x-element columns across y and all physical Fourier
    // planes. Full-span z integration preserves the spectral SVV work identity.
    for (int plane = 0; plane < planes; ++plane)
    {
        auto planeGeometry = field->GetPlane(plane);
        for (int element = 0; element < planeGeometry->GetExpSize(); ++element)
        {
            auto expansion = planeGeometry->GetExp(element);
            const int count = expansion->GetTotPoints();
            const int offset = plane * planePoints + planeGeometry->GetPhys_Offset(element);
            Array<OneD, NekDouble> localOnes(count, 1.0), energy(count, 0.0),
                divergenceSquared(count);
            for (int q = 0; q < count; ++q)
            {
                const int point = offset + q;
                for (int component = 0; component < 3; ++component)
                {
                    const auto velocity = physical[component][point];
                    const auto fluctuation = velocity - (component == 0 ? mean : 0.0);
                    energy[q] += 0.5 * fluctuation * fluctuation;
                }
                divergenceSquared[q] = divergence[point] * divergence[point];
            }
            auto &integrals = slabs.at(intervals[element]);
            integrals[0] += weight * expansion->Integral(localOnes);
            for (int component = 0; component < 3; ++component)
            {
                integrals[component + 1] += weight * expansion->Integral(physical[component] + offset);
            }
            integrals[4] += weight * expansion->Integral(energy);
            integrals[5] += 2 * viscosity * weight * expansion->Integral(strainSquared + offset);
            integrals[6] += weight * expansion->Integral(svvDensity[0] + offset);
            integrals[7] += weight * expansion->Integral(svvDensity[1] + offset);
            integrals[8] += weight * expansion->Integral(divergenceSquared);
        }
    }
    std::ofstream output(session->GetSessionName() + ".slabs.csv");
    output << std::setprecision(17)
           << "xmin,xmax,xmid,volume,mean_u,mean_v,mean_w,"
           << "instantaneous_K_about_prescribed_mean,molecular_strain_dissipation,"
           << "svv_xy,svv_z,divergence_rms\n";
    NekDouble totalVolume = 0.0;
    for (const auto &slab : slabs)
    {
        const auto &integrals = slab.second;
        const auto slabVolume = integrals[0];
        const auto expectedVolume = (slab.first.second - slab.first.first) *
                                    (yMax - yMin) * field->GetHomoLen();
        ASSERTL0(slabVolume > 0 && std::abs(slabVolume - expectedVolume) <=
                     1e-9 * expectedVolume,
                 "Slab profiles require complete y columns and full-span z volume");
        for (const auto value : integrals)
        {
            ASSERTL0(std::isfinite(value), "Non-finite slab integral");
        }
        output << slab.first.first << ',' << slab.first.second << ','
               << 0.5 * (slab.first.first + slab.first.second) << ',' << slabVolume;
        // K is instantaneous energy about fixed Ubar, not temporal-central TKE.
        for (int quantity = 1; quantity < 8; ++quantity)
        {
            output << ',' << integrals[quantity] / slabVolume;
        }
        output << ',' << std::sqrt(std::max(0.0, integrals[8] / slabVolume)) << '\n';
        totalVolume += slabVolume;
    }
    ASSERTL0(std::abs(totalVolume - volume) <= 1e-10 * volume,
             "Slab volumes do not cover the imported field");
    ASSERTL0(output.good(), "Cannot write HIT slab profiles");
}
} // namespace

int main(int argc, char *argv[])
{
    if (argc < 3)
    {
        std::cerr << "Usage: HITFlowDiagnostics [session options] session.xml solution.fld\n";
        return 2;
    }
    try
    {
        std::vector<char *> arguments(argv, argv + argc - 1);
        arguments.push_back(nullptr);
        auto session = LibUtilities::SessionReader::CreateInstance(
            argc - 1, arguments.data());
        auto graph = SpatialDomains::MeshGraphIO::Read(session);
        auto fields = HITTestTools::Fields(session, graph);
        auto data = HITTestTools::ImportFields(session, fields, argv[argc - 1]);
        const int points = fields[0]->GetNpoints();
        const bool homogeneous = fields[0]->GetExpType() == MultiRegions::e3DH1D;
        Array<OneD, NekDouble> divergence(points, 0.0), ones(points, 1.0);
        std::array<std::array<Array<OneD, NekDouble>, 3>, 3> derivatives;
        NekDouble velocitySquared = 0.0, gradientSquared = 0.0;
        for (int component = 0; component < 3; ++component)
        {
            Array<OneD, NekDouble> dx(points), dy(points), dz(points);
            fields[component]->PhysDeriv(data.physical[component], dx, dy, dz);
            derivatives[component] = {dx, dy, dz};
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
        NekDouble viscosity;
        session->LoadParameter("Kinvis", viscosity);
        Array<OneD, NekDouble> strainSquared(points, 0.0);
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                for (int q = 0; q < points; ++q)
                {
                    const NekDouble strain = 0.5 *
                        (derivatives[i][j][q] + derivatives[j][i][q]);
                    strainSquared[q] += strain * strain;
                }
            }
        }
        // Instantaneous total resolved strain, not a temporal fluctuation
        // statistic. Its distinction from Reynolds-stress dissipation matters
        // during startup and in a spatially developing mean field.
        const NekDouble resolvedStrainDissipation =
            2.0 * viscosity * fields[0]->Integral(strainSquared) / volume;
        std::array<Array<OneD, NekDouble>, 2> svvDensity;
        const auto svvWork = HITTestTools::SVVWork(session, fields[0], derivatives,
            volume, homogeneous ? &svvDensity : nullptr);
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
        if (homogeneous)
        {
            WriteSlabs(session, fields[0], data.physical, strainSquared,
                       divergence, svvDensity, viscosity, volume);
        }
        std::ofstream output(session->GetSessionName() + ".flow.json");
        output << std::setprecision(17)
               << "{\"volume\":" << volume
               << ",\"velocity_l2\":" << velocityNorm
               << ",\"reference_length\":" << referenceLength
               << ",\"gradient_l2\":" << std::sqrt(gradientSquared)
               << ",\"viscosity\":" << viscosity
               << ",\"instantaneous_resolved_strain_dissipation\":"
               << resolvedStrainDissipation
               << ",\"instantaneous_svv_work_xy\":" << svvWork[0]
               << ",\"instantaneous_svv_work_z\":" << svvWork[1]
               << ",\"divergence_rms\":" << divNorm / std::sqrt(volume)
               << ",\"divergence_fixed_cbc_reference\":"
               << 0.024 * divNorm / (0.222 * std::sqrt(volume))
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
