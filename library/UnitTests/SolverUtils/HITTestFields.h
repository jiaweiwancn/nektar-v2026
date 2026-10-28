///////////////////////////////////////////////////////////////////////////////
// File: HITTestFields.h
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_UNITTESTS_HITTESTFIELDS_H
#define NEKTAR_UNITTESTS_HITTESTFIELDS_H

#include <MultiRegions/ContField.h>
#include <MultiRegions/ContField3DHomogeneous1D.h>
#include <SpatialDomains/MeshGraphIO.h>

namespace Nektar::HITTestTools
{
inline Array<OneD, MultiRegions::ExpListSharedPtr> Fields(
    const LibUtilities::SessionReaderSharedPtr &session,
    const SpatialDomains::MeshGraphSharedPtr &graph)
{
    Array<OneD, MultiRegions::ExpListSharedPtr> result(4);
    const bool homogeneous = session->DefinesSolverInfo("HOMOGENEOUS");
    bool useFFT = false;
    session->MatchSolverInfo("USEFFT", "FFTW", useFFT, false);
    for (int variable = 0; variable < 4; ++variable)
    {
        if (homogeneous)
        {
            const int modes = session->GetParameter("HomModesZ");
            const LibUtilities::PointsKey points(modes,
                LibUtilities::eFourierEvenlySpaced);
            const LibUtilities::BasisKey basis(LibUtilities::eFourier,
                                               modes, points);
            result[variable] =
                MemoryManager<MultiRegions::ContField3DHomogeneous1D>::
                    AllocateSharedPtr(session, basis,
                        session->GetParameter("LZ"), useFFT, false, graph,
                        session->GetVariable(variable), true);
        }
        else
        {
            result[variable] = MemoryManager<MultiRegions::ContField>::
                AllocateSharedPtr(session, graph,
                    session->GetVariable(variable), true, true);
        }
    }
    return result;
}
} // namespace Nektar::HITTestTools

#endif
