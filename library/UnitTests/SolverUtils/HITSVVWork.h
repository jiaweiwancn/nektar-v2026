///////////////////////////////////////////////////////////////////////////////
// File: HITSVVWork.h
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_HITSVVWORK_H
#define NEKTAR_HITSVVWORK_H

#include <MultiRegions/ExpList.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <regex>

namespace Nektar::HITTestTools
{
// Instantaneous work of the configured positive SVV operator. This diagnostic
// is restricted to the affine H1D pilot boxes with constant prescribed SVV speed.
inline std::array<NekDouble, 2> SVVWork(
    const LibUtilities::SessionReaderSharedPtr &session,
    const MultiRegions::ExpListSharedPtr &field,
    const std::array<std::array<Array<OneD, NekDouble>, 3>, 3> &derivatives,
    NekDouble volume)
{
    std::array<NekDouble, 2> work{};
    ASSERTL0(std::isfinite(volume) && volume > 0.0,
             "SVV work requires positive finite volume");
    bool dg = false, exponential = false, fourier = false;
    session->MatchSolverInfo("SpectralVanishingViscositySpectralHP", "DGKernel", dg, false);
    session->MatchSolverInfo("SpectralVanishingViscositySpectralHP", "ExpKernel", exponential, false);
    session->MatchSolverInfo("SpectralVanishingViscosityHomo1D", "ExpKernel", fourier, false);
    ASSERTL0(!session->DefinesSolverInfo("SpectralVanishingViscosity") &&
                 (!session->DefinesSolverInfo("SpectralVanishingViscositySpectralHP") || dg || exponential) &&
                 (!session->DefinesSolverInfo("SpectralVanishingViscosityHomo1D") || fourier),
             "Pilot SVV work supports explicit SpectralHP DG/Exp and Homo1D Exp switches only");
    if (!dg && !exponential && !fourier)
    {
        return work;
    }
    ASSERTL0(field->GetExpType() == MultiRegions::e3DH1D &&
                 session->GetComm()->GetSize() == 1,
             "Pilot SVV work requires a complete serial H1D field");
    const auto ids = field->GetZIDs();
    const int planes = ids.size();
    const int planePoints = field->GetPlane(0)->GetNpoints();
    const auto length = field->GetHomoLen();
    NekDouble coefficient, cutoff;
    session->LoadParameter("SVVDiffCoeff", coefficient, dg ? 1.0 : 0.1);
    session->LoadParameter("SVVCutoffRatio", cutoff, 0.75);
    if (dg || exponential)
    {
        NekDouble speed = 1.0;
        if (dg && session->DefinesFunction("SVVVelocityMagnitude"))
        {
            speed = 0.0;
            for (const auto &variable : {"u", "v", "w"})
            {
                const auto function = session->GetFunction("SVVVelocityMagnitude", variable);
                ASSERTL0(!std::regex_search(function->GetExpression(),
                             std::regex("\\b(x|y|z|t|awgn|rand|random)\\b", std::regex::icase)),
                         "Pilot SVV work requires constant SVVVelocityMagnitude");
                const auto component = function->Evaluate();
                speed += component * component;
            }
            speed = std::sqrt(speed);
        }
        for (int plane = 0; plane < planes; ++plane)
        {
            auto geometry = field->GetPlane(plane);
            for (int element = 0; element < geometry->GetExpSize(); ++element)
            {
                auto expansion = geometry->GetExp(element);
                ASSERTL0(expansion->DetShapeType() == LibUtilities::eQuadrilateral &&
                             expansion->GetGeomFactors()->GetGtype() == SpatialDomains::eRegular,
                         "Pilot SVV work supports affine quadrilaterals");
                NekDouble h = 0.0;
                for (int edge = 0; edge < expansion->GetGeom()->GetNumEdges(); ++edge)
                {
                    const auto segment = expansion->GetGeom()->GetEdge(edge);
                    h = std::max(h, segment->GetVertex(0)->dist(*segment->GetVertex(1)));
                }
                const int p = std::max(expansion->GetBasisNumModes(0),
                                       expansion->GetBasisNumModes(1)) - 1;
                StdRegions::ConstFactorMap factors;
                factors[StdRegions::eFactorSVVCutoffRatio] = cutoff;
                factors[StdRegions::eFactorSVVDiffCoeff] = coefficient;
                if (dg)
                {
                    factors[StdRegions::eFactorSVVDGKerDiffCoeff] = speed * h / p;
                }
                const StdRegions::StdMatrixKey key(StdRegions::eLaplacian,
                    expansion->DetShapeType(), *expansion, factors);
                const int count = expansion->GetTotPoints();
                const int offset = plane * planePoints + geometry->GetPhys_Offset(element);
                for (int component = 0; component < 3; ++component)
                {
                    for (int direction = 0; direction < 2; ++direction)
                    {
                        Array<OneD, NekDouble> filtered(count), product(count);
                        for (int q = 0; q < count; ++q)
                        {
                            filtered[q] = derivatives[component][direction][offset + q];
                        }
                        expansion->SVVLaplacianFilter(filtered, key);
                        for (int q = 0; q < count; ++q)
                        {
                            product[q] = filtered[q] * derivatives[component][direction][offset + q];
                        }
                        work[0] += length / planes * expansion->Integral(product) / volume;
                    }
                }
            }
        }
    }
    if (fourier)
    {
        NekDouble zCoefficient, zCutoff;
        session->LoadParameter("SVVDiffCoeffHomo1D", zCoefficient, coefficient);
        session->LoadParameter("SVVCutoffRatioHomo1D", zCutoff, cutoff);
        const int maximum = planes / 2;
        const int start = zCutoff * maximum;
        for (int component = 0; component < 3; ++component)
        {
            Array<OneD, NekDouble> transformed(field->GetNpoints());
            field->HomogeneousFwdTrans(field->GetNpoints(), derivatives[component][2], transformed);
            NekDouble parseval = 0.0;
            for (int plane = 0; plane < planes; ++plane)
            {
                Array<OneD, NekDouble> squared(planePoints);
                for (int q = 0; q < planePoints; ++q)
                {
                    const auto value = transformed[plane * planePoints + q];
                    squared[q] = value * value;
                }
                const auto norm = field->GetPlane(plane)->Integral(squared) *
                                  (ids[plane] == 0 ? length : length / 2.0);
                parseval += norm;
                const int wave = ids[plane] / 2;
                if (wave > start)
                {
                    const NekDouble ratio = NekDouble(wave - maximum) / (wave - start);
                    work[1] += zCoefficient * std::exp(-ratio * ratio) * norm / volume;
                }
            }
            const auto physicalNorm = field->L2(derivatives[component][2]);
            ASSERTL0(std::abs(parseval - physicalNorm * physicalNorm) <=
                         1e-10 * std::max(parseval, physicalNorm * physicalNorm),
                     "Fourier SVV diagnostic failed its Parseval check");
        }
    }
    ASSERTL0(std::isfinite(work[0]) && std::isfinite(work[1]) &&
                 work[0] >= 0.0 && work[1] >= 0.0,
             "Non-finite or negative SVV work");
    return work;
}
} // namespace Nektar::HITTestTools
#endif
