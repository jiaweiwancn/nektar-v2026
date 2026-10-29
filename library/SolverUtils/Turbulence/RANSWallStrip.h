///////////////////////////////////////////////////////////////////////////////
// File: RANSWallStrip.h
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_RANS_WALL_STRIP_H
#define NEKTAR_SOLVERUTILS_RANS_WALL_STRIP_H

#include <cmath>
#include <stdexcept>

namespace Nektar::SolverUtils
{
/**
 * The scalar, per-unit-wall-area contract for a flat control-volume wall strip.
 * The CG domain starts at H, not at the physical wall. P is a strip state and
 * N is a physical sample in the CG domain. These locations are independent of
 * quadrature and polynomial order. This is an assembly primitive, not a RANS
 * closure, a wall function, or an installed incompressible solver option.
 */
struct RANSWallStripGeometry
{
    double height, pointP, pointN;

    void Validate() const
    {
        if (!(std::isfinite(height) && std::isfinite(pointP) &&
              std::isfinite(pointN) && 0.0 < pointP && pointP < height &&
              height < pointN))
        {
            throw std::invalid_argument("RANS wall strip requires 0 < P < H < N");
        }
    }

    // Log-coordinate interpolation reproduces C + a log(z) exactly. It avoids
    // equating the physical CG trace at H to an inconsistent linear-z estimate.
    double VelocityWeightN() const
    {
        Validate();
        return std::log(height / pointP) / std::log(pointN / pointP);
    }

    double KWeightN() const
    {
        Validate();
        return (height - pointP) / (pointN - pointP);
    }

    double VelocityAtH(double valueP, double valueN) const
    {
        const double weight = VelocityWeightN();
        return (1.0 - weight) * valueP + weight * valueN;
    }

    double KAtH(double valueP, double valueN) const
    {
        const double weight = KWeightN();
        return (1.0 - weight) * valueP + weight * valueN;
    }

    double EliminateVelocityP(double valueH, double valueN) const
    {
        const double weight = VelocityWeightN();
        return (valueH - weight * valueN) / (1.0 - weight);
    }

    double EliminateKP(double valueH, double valueN) const
    {
        const double weight = KWeightN();
        return (valueH - weight * valueN) / (1.0 - weight);
    }
};

struct RANSWallStripState
{
    double velocity, k, dissipation;
};

struct RANSWallStripRates
{
    double velocity, k;
};

struct RANSWallStripTrace
{
    double velocityH, velocityN, kH, kN, dissipationH;
};

struct RANSWallStripExchange
{
    // Positive z is away from the wall. These are gradient fluxes inside the
    // diffusion operator. The CG outward flux at its bottom is their negative.
    double momentumH, tractionWall, kH, kWall;
    double meanMomentumForcing, meanKAdvection;
    double meanProduction, meanDissipation;
    // Both values come from the SAME selected auxiliary dissipation profile.
    double prescribedDissipationP, prescribedDissipationH;
};

struct RANSWallStripResidual
{
    double momentum, k, dissipationP, velocityMatch, kMatch, dissipationH;
};

inline RANSWallStripResidual AssembleRANSWallStrip(
    const RANSWallStripGeometry &geometry, const RANSWallStripState &state,
    const RANSWallStripRates &rate, const RANSWallStripTrace &trace,
    const RANSWallStripExchange &exchange)
{
    geometry.Validate();
    const double H = geometry.height;
    return {
        H * (rate.velocity + exchange.meanMomentumForcing) -
            exchange.momentumH + exchange.tractionWall,
        H * (rate.k + exchange.meanKAdvection) - exchange.kH + exchange.kWall -
            H * (exchange.meanProduction - exchange.meanDissipation),
        state.dissipation - exchange.prescribedDissipationP,
        trace.velocityH - geometry.VelocityAtH(state.velocity, trace.velocityN),
        trace.kH - geometry.KAtH(state.k, trace.kN),
        trace.dissipationH - exchange.prescribedDissipationH};
}

// Integrated incompressibility for a flat strip. In multiple dimensions the
// divergence is that of the strip's tangential velocity, in its physical frame.
inline double RANSWallStripMassResidual(
    const RANSWallStripGeometry &geometry, double normalVelocityH,
    double normalVelocityWall, double tangentialDivergence)
{
    geometry.Validate();
    return normalVelocityH - normalVelocityWall +
           geometry.height * tangentialDivergence;
}
} // namespace Nektar::SolverUtils

#endif
