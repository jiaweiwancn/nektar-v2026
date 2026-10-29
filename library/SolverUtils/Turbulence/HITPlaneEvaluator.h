///////////////////////////////////////////////////////////////////////////////
// File: HITPlaneEvaluator.h
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_HITPLANEEVALUATOR_H
#define NEKTAR_SOLVERUTILS_HITPLANEEVALUATOR_H

#include <SolverUtils/Turbulence/HITGenerator.h>
#include <string>

namespace Nektar::SolverUtils
{
/// Optional velocity-only acceleration at immutable points on one x plane.
/// The caller must keep the generator alive. Sampling remains an absolute-time
/// function; there is no phase recurrence, random draw, or mutable time state.
class SOLVER_UTILS_EXPORT HITPlaneEvaluator
{
public:
    HITPlaneEvaluator(const HITGenerator &generator,
                      const std::vector<HITVector> &points,
                      std::size_t maxCacheBytes = 64 * 1024 * 1024);

    void Evaluate(NekDouble time, std::vector<HITVector> &velocities) const;

    /// Additional vector payload, including per-call accumulation workspace.
    /// Excludes allocator overhead, constructor maps, input points and generator.
    std::size_t CacheBytes() const;
    bool IsGrouped() const;
    const std::string &Method() const;

private:
    struct Term
    {
        const HITFourierMode *mode;
        std::size_t xPhase;
        std::size_t positiveBin;
        std::size_t negativeBin;
    };

    bool SafeAtTime(NekDouble time) const;
    void Reference(NekDouble time, std::vector<HITVector> &velocities) const;

    const HITGenerator &m_generator;
    std::vector<HITVector> m_points;
    HITVector m_maxAbsWave{};
    std::vector<Term> m_terms;
    std::vector<NekDouble> m_xWaves;
    std::vector<std::array<NekDouble, 2>> m_transverseWaves;
    std::vector<std::complex<NekDouble>> m_spatialPhases;
    std::size_t m_cacheBytes = 0;
    bool m_grouped = false;
    std::string m_method = "Reference";
};
} // namespace Nektar::SolverUtils

#endif
