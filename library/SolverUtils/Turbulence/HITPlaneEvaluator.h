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
    /// Excludes allocator overhead, input points and generator.
    std::size_t CacheBytes() const;
    bool IsGrouped() const;
    const std::string &Method() const;

private:
    bool SafeAtTime(NekDouble time) const;
    std::size_t TransverseBin(int y, int z) const;
    void EvaluateSeparable(const std::vector<HITComplexVector> &amplitude,
                           std::vector<HITVector> &velocities) const;
    void Reference(NekDouble time, std::vector<HITVector> &velocities) const;

    const HITGenerator &m_generator;
    std::vector<HITVector> m_points;
    HITVector m_maxAbsWave{};
    // Mode indices recover their phase/bin positions from these bounded
    // lattice tables, without an auxiliary record for every retained 3D mode.
    std::vector<std::size_t> m_xIndex;
    std::vector<std::size_t> m_transverseIndex;
    std::array<std::size_t, 2> m_transverseOffset{};
    std::size_t m_transverseStride = 0;
    std::vector<NekDouble> m_xWaves;
    std::vector<std::array<NekDouble, 2>> m_transverseWaves;
    std::vector<std::complex<NekDouble>> m_spatialPhases;
    std::array<std::vector<std::complex<NekDouble>>, 2> m_axisPhases;
    std::vector<std::array<std::size_t, 2>> m_pointIndex;
    std::array<std::size_t, 2> m_uniqueCoordinates{};
    std::size_t m_firstAxis = 0;
    std::size_t m_cacheBytes = 0;
    bool m_separable = false;
    bool m_grouped = false;
    std::string m_method = "Reference";
};
} // namespace Nektar::SolverUtils

#endif
