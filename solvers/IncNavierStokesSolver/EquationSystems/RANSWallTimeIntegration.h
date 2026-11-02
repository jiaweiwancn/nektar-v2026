///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_RANS_WALL_TIME_INTEGRATION_H
#define NEKTAR_RANS_WALL_TIME_INTEGRATION_H

#include <IncNavierStokesSolver/EquationSystems/RANSHistory.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallFieldLayout.h>
#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace Nektar
{
/** Native constant-step IMEX for a flattened differential coefficient vector.
 *
 * A single native variable accommodates the different native coefficient
 * counts of velocity, full-domain k and outer-only epsilon. Algebraic pressure
 * and gauge values stay in the caller's coupled stage state. The implicit
 * callback solves the complete coupled stage for history and lambda, with
 * alpha=1/lambda, and returns its differential coefficients. Its input already
 * contains the native explicit-history combination; it is not a previous-time
 * solution to be reconstructed by the caller.
 *
 * The original RANSIMEX tableaux, order-1 startup, explicit derivative slots
 * and GLM output combinations are retained. All algebraic/essential/dependent
 * constraints are the coupled callback's responsibility. The native projection
 * callback is the identity, including native startup phase conversion.
 */
class RANSWallTimeIntegration
{
public:
    using Vector = RANSWallAlgebra::Vector;
    using ImplicitCallback =
        std::function<Vector(const Vector &, double, double)>;
    using ExplicitCallback = std::function<Vector(const Vector &, double)>;

    struct Snapshot
    {
        std::size_t order = 0, steps = 0;
        double time = 0, dt = 0;
        Vector state;
        // Native order-N IMEX stores N solution values, followed by N
        // dt-scaled explicit derivatives, in this exact native order.
        std::vector<Vector> history;
        Vector times;
    };

    struct Checkpoint
    {
        std::string signature;
        Snapshot integration;
        Vector coupledState;
        Vector fixedReference;
        Vector fixedDefect;
        double eta = 0;
    };

    RANSWallTimeIntegration(std::size_t order, double dt,
                           ImplicitCallback implicit,
                           ExplicitCallback explicitRhs);
    RANSWallTimeIntegration(const RANSWallTimeIntegration &) = delete;
    RANSWallTimeIntegration &operator=(const RANSWallTimeIntegration &) = delete;
    RANSWallTimeIntegration(RANSWallTimeIntegration &&) = delete;
    RANSWallTimeIntegration &operator=(RANSWallTimeIntegration &&) = delete;

    void Initialise(const Vector &state, double time);
    Vector Advance();
    void Restore(const Snapshot &snapshot);
    Snapshot Capture() const;

    const Vector &DifferentialState() const { return m_state; }
    double Time() const { return m_time; }
    double TimeStep() const { return m_dt; }
    std::size_t Steps() const { return m_steps; }
    std::size_t Order() const { return m_order; }

    /** Publish a new checkpoint directory without overwriting any generation.
     * Only the caller's selected writer rank invokes this method. The returned
     * directory contains one atomically renamed State.bin; incomplete writes
     * cannot be read as complete checkpoints. All vectors use the caller's
     * canonical global ordering and are independent of MPI partition count.
     */
    static std::string WriteCheckpoint(const std::string &baseDirectory,
                                       const Checkpoint &checkpoint);
    static Checkpoint ReadCheckpoint(const std::string &directory,
                                     const std::string &expectedSignature,
                                     std::size_t expectedOrder, double expectedDt,
                                     std::size_t expectedCoupledSize,
                                     std::size_t expectedDifferentialSize);

private:
    using Operators = LibUtilities::TimeIntegrationSchemeOperators;
    void BuildScheme(const Vector &initial, double time,
                     const Snapshot *restored);
    void Implicit(Operators::InArrayType &input, Operators::OutArrayType &output,
                  double time, double lambda);
    void Explicit(Operators::InArrayType &input, Operators::OutArrayType &output,
                  double time);
    void Projection(Operators::InArrayType &input,
                    Operators::OutArrayType &output, double time);
    Vector Input(Operators::InArrayType &input) const;
    void Output(const Vector &value, Operators::OutArrayType &output) const;

    std::size_t m_order, m_steps = 0, m_localSteps = 0, m_size = 0;
    double m_dt, m_time = 0;
    ImplicitCallback m_implicit;
    ExplicitCallback m_explicit;
    std::shared_ptr<RANSIMEX> m_scheme;
    Vector m_state;
};
}
#endif
