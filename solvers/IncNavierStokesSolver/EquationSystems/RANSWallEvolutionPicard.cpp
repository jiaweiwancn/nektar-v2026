///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric>

namespace Nektar
{
using namespace SolverUtils;
using RANSWallAlgebra::Maximum;
using RANSWallAlgebra::Multiply;
using RANSWallAlgebra::Require;
namespace J = RANSWallSolveJson;

RANSWallEvolutionPatch::PicardBlock RANSWallEvolutionPatch::AssemblePicardBlock(
    int block, const Vector &s, const Vector &history, double alpha,
    double time, const RANSWallEvolutionData &data,
    const std::map<int, double> &external) const
{
    const bool flow = block == 0;
    PicardBlock system;
    std::vector<int> index(Size(), -1);
    for (int row = 0; row < Size(); ++row)
    {
        const bool selected =
            flow ? row < offset[K()] || row >= DifferentialSize()
                 : row >= offset[block] && row < offset[block + 1];
        if (selected && !external.count(row))
        {
            index[row] = system.rows.size();
            system.rows.push_back(row);
        }
    }
    const auto n = system.rows.size();
    system.matrix.assign(n, Vector(n));
    system.residual.assign(n, 0.);
    const auto &model = m_parameters.model.Parameters();
    const double nu   = m_parameters.geometry.molecularNu;
    Refresh(s, time, &data);

    auto addScalar = [&](const RANSWallGlobalPoint &point,
                         const RANSVector &velocity, double weight,
                         double diffusion, double rate, double source) {
        const double value = point.Value(s, offset[block]);
        const auto grad    = point.Gradient(s, offset[block]);
        for (int i = 0; i < point.ids.size(); ++i)
        {
            const int row = index[offset[block] + point.ids[i]];
            if (row < 0)
                continue;
            const double test = point.phi[i];
            const double adv  = point.derivative[0][i] * velocity[0] +
                               point.derivative[2][i] * velocity[2];
            const double dx = point.derivative[0][i],
                         dz = point.derivative[2][i];
            system.residual[row] +=
                weight * (test * (rate * value - source) - adv * value +
                          diffusion * (dx * grad[0] + dz * grad[2]));
            for (int j = 0; j < point.ids.size(); ++j)
            {
                const int col = index[offset[block] + point.ids[j]];
                if (col >= 0)
                    system.matrix[row][col] +=
                        weight * ((test * rate - adv) * point.phi[j] +
                                  diffusion * (dx * point.derivative[0][j] +
                                               dz * point.derivative[2][j]));
            }
        }
    };
    auto evaluate = [&](const VolumePoint &q, const Ray *ray,
                        const RANSVector &force,
                        const RANSWallGlobalPoint *forcingPressure) {
        const bool inner = ray != nullptr;
        RANSVector u{};
        RANSTensor grad{};
        for (int a = 0; a < m_dimension; ++a)
        {
            const int axis = PhysicalAxis(a);
            u[axis]        = q.v.Value(s, offset[a]);
            grad[axis]     = q.v.Gradient(s, offset[a]);
        }
        const double k = q.v.Value(s, offset[K()]);
        double viscosity;
        RANSTensor stress{};
        const int pBlock = inner ? InnerPressure() : OuterPressure();
        if (inner)
        {
            const double nt = ray->profile->Original().EddyViscosity(
                Distance(q.v.coordinate));
            const auto law = InnerStress(k, nt, grad);
            if (!flow)
            {
                double rate = 0.;
                WallSource(s, q.v.coordinate, law.production, *ray->profile,
                           &ray->destruction, &rate);
                addScalar(q.v, u, q.weight, nu + nt / model.sigmaK, rate,
                          law.production +
                              data.innerKSource(q.v.coordinate, time));
                return;
            }
            viscosity = nu + law.tangentialTangent;
            stress =
                RANSWallInnerMeanStress(law, nu, q.p.Value(s, offset[pBlock]));
        }
        else
        {
            const double epsilon = q.e.Value(s, offset[Epsilon()]);
            const auto law = OuterClosure(k, epsilon, grad, q.v.coordinate);
            if (!flow)
            {
                if (block == K())
                {
                    const double rate =
                        (epsilon - std::min(0., law.production)) / k;
                    addScalar(q.v, u, q.weight, law.kDiffusion, rate,
                              std::max(0., law.production) +
                                  data.outerKSource(q.v.coordinate, time));
                }
                else
                {
                    addScalar(
                        q.e, u, q.weight, law.epsilonDiffusion,
                        model.cEpsilon2 * epsilon / k,
                        law.epsilonProduction + law.yap +
                            data.outerEpsilonSource(q.v.coordinate, time));
                }
                return;
            }
            viscosity = law.eddyViscosity + m_regularizationGamma +
                        m_regularizationScale * m_parameters.model.Kappa() *
                            m_robustUstar * Distance(q.v.coordinate);
            const double pi = q.p.Value(s, offset[pBlock]);
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                    stress[a][b] = -law.reynoldsStress[a][b] +
                                   (a == b ? 2 * k / 3 - pi : 0.);
        }
        const double divergence = grad[0][0] + grad[2][2];
        for (int p = 0; p < q.p.ids.size(); ++p)
            system.residual[index[offset[pBlock] + q.p.ids[p]]] +=
                q.weight * q.p.phi[p] * divergence;

        // Oseen advection and a positive symmetric-strain part are implicit.
        // The load at s includes the complete stress, so its remainder
        // (including CLS nonlinear stresses and the reference-strain load)
        // is deferred. Both regions use the shared velocity test space.
        for (int a = 0; a < m_dimension; ++a)
        {
            const int aa = PhysicalAxis(a);
            for (int i = 0; i < q.v.ids.size(); ++i)
            {
                const int row               = index[offset[a] + q.v.ids[i]];
                const double dx             = q.v.derivative[0][i],
                             dz             = q.v.derivative[2][i];
                const double testDerivative = q.v.derivative[aa][i];
                if (row >= 0)
                {
                    system.residual[row] +=
                        q.weight * (dx * (stress[aa][0] - u[aa] * u[0]) +
                                    dz * (stress[aa][2] - u[aa] * u[2]) -
                                    q.v.phi[i] * force[aa]);
                    for (int b = 0; b < m_dimension; ++b)
                    {
                        const int bb = PhysicalAxis(b);
                        for (int j = 0; j < q.v.ids.size(); ++j)
                        {
                            const int col = index[offset[b] + q.v.ids[j]];
                            if (col < 0)
                                continue;
                            double value = viscosity * q.v.derivative[bb][i] *
                                           q.v.derivative[aa][j];
                            if (a == b)
                                value +=
                                    viscosity * (dx * q.v.derivative[0][j] +
                                                 dz * q.v.derivative[2][j]) -
                                    (dx * u[0] + dz * u[2]) * q.v.phi[j];
                            system.matrix[row][col] += q.weight * value;
                        }
                    }
                    if (forcingPressure && aa == 0)
                    {
                        // Flat tangential pressure gradient in the roughness
                        // load. The mechanical -2 grad(k)/3 and convection
                        // parts of CartesianForcing stay frozen in this solve.
                        const auto &fp = *forcingPressure;
                        for (int p = 0; p < fp.ids.size(); ++p)
                            system.matrix[row][index[offset[OuterPressure()] +
                                                     fp.ids[p]]] -=
                                q.weight * q.v.phi[i] * fp.derivative[0][p];
                    }
                }
                if (row < 0)
                    continue;
                for (int p = 0; p < q.p.ids.size(); ++p)
                {
                    const int pressure = index[offset[pBlock] + q.p.ids[p]];
                    const double entry = q.weight * q.p.phi[p] * testDerivative;
                    system.matrix[row][pressure] -= entry;
                    system.matrix[pressure][row] += entry;
                }
            }
        }
    };

    if (block != Epsilon())
        for (const auto &ray : m_rays)
        {
            RANSVector force{};
            if (flow)
                force = CartesianForcing(s, ray.coordinate, time, data, false);
            for (const auto &q : ray.points)
            {
                const bool rough = flow && Distance(q.v.coordinate) <=
                                               m_parameters.geometry.roughness;
                evaluate(q, &ray, rough ? force : RANSVector{},
                         rough ? &ray.outerP : nullptr);
            }
        }
    for (const auto &q : m_outer)
        evaluate(q, nullptr, {}, nullptr);

    const auto body = flow ? OwnedBodyLoad(time, data) : Vector{};
    for (int b = flow ? 0 : block; b < (flow ? K() : block + 1); ++b)
    {
        Vector change(size[b]);
        for (int j = 0; j < size[b]; ++j)
            change[j] = s[offset[b] + j] - history[offset[b] + j];
        const auto storage = Multiply(m_mass[b], change);
        for (int i = 0; i < size[b]; ++i)
        {
            const int row = index[offset[b] + i];
            if (row < 0)
                continue;
            system.residual[row] += alpha * storage[i];
            if (flow)
                system.residual[row] -= body[offset[b] + i];
            if (m_eta)
                system.residual[row] -= m_eta * m_defect[offset[b] + i];
            for (int j = 0; j < size[b]; ++j)
            {
                const int col = index[offset[b] + j];
                if (col >= 0)
                    system.matrix[row][col] +=
                        (alpha + m_picardPseudoAlpha) * m_mass[b][i][j];
            }
        }
    }
    if (flow)
    {
        const int gauge = index[offset[Gauge()]];
        for (int side = 0; side < 2; ++side)
            for (int i = 0; i < m_pressureMean[side].size(); ++i)
            {
                const int id              = offset[InnerPressure() + side] + i,
                          row             = index[id];
                const double w            = m_pressureMean[side][i];
                system.matrix[row][gauge] = system.matrix[gauge][row] = w;
                system.residual[row] += s[offset[Gauge()]] * w;
                system.residual[gauge] += s[id] * w;
            }
    }
    else if (block == Epsilon())
    {
        // Hold Phi fixed for this scalar solve. Free transport rows retain
        // the columns of the matching trace, including its update.
        const auto trace = DissipationTrace(s, time, data);
        for (const auto &[id, value] : trace[0])
        {
            const int row = index[offset[block] + id];
            Require(row >= 0,
                    "Picard matching trace overlaps external lifting");
            std::fill(system.matrix[row].begin(), system.matrix[row].end(), 0.);
            system.matrix[row][row] = 1.;
            system.residual[row] =
                s[offset[block] + id] - (value + trace[1].at(id));
        }
    }
    return system;
}

RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::SolveSegregatedPicard(
    Vector s, const Vector &history, double alpha, double time,
    const RANSWallEvolutionData &data, int &iterations, double tolerance,
    int maximumIterations) const
{
    Require(s.size() == Size() && history.size() == s.size() &&
                std::isfinite(alpha) && alpha > 0. && std::isfinite(time) &&
                std::isfinite(tolerance) && tolerance > 0. &&
                maximumIterations >= 0,
            "invalid Picard physical stage");
    iterations = 0;
    m_lastFeasibleState.clear();
    const auto incomingGuessHash = RANSWallR3::Hash(s);
    s                            = ExternalLift(s, time, data);
    auto external                = ExternalConstraints(time, data);
    for (int id : m_traceIds)
        external.erase(offset[Epsilon()] + id);
    Require(external.empty() || external.rbegin()->first < DifferentialSize(),
            "Picard external lifting may constrain differential rows only");
    const auto start = std::chrono::steady_clock::now();
    int flowSolves = 0, kSolves = 0, epsilonSolves = 0;
    auto emit = [&](const std::string &event,
                    std::vector<std::pair<std::string, std::string>> values) {
        values.emplace_back("event", J::Quote(event));
        values.emplace_back("sweep", J::Number(iterations));
        values.emplace_back("alpha", J::Number(alpha));
        values.emplace_back("time", J::Number(time));
        values.emplace_back(
            "seconds", J::Number(std::chrono::duration<double>(
                                     std::chrono::steady_clock::now() - start)
                                     .count()));
        std::cout << "EP-AWF Picard " << J::Object(values) << std::endl;
    };
    auto admissible = [&](const Vector &state, ScalarDomainResult &domain,
                          std::string &ground) {
        domain     = ScalarDomain(state, time, data);
        bool valid = false;
        if (domain.admissible)
            ground = WallTraceJSON(domain.effectiveState, valid);
        return valid;
    };
    auto normsAt = [&](const Vector &state) {
        const auto residual = Residual(state, history, alpha, time, data);
        return Norms(residual, Constraints(state, time, data), m_robustUstar);
    };
    auto trialNorms = [&](const Vector &state, double fraction) {
        try
        {
            return normsAt(state);
        }
        catch (const std::exception &error)
        {
            emit("rejected_residual", {{"fraction", J::Number(fraction)},
                                       {"reason", J::Quote(error.what())}});
            return Vector{};
        }
    };
    auto blockCheck = [&](const PicardBlock &system, const Vector &state,
                          const char *name) {
        const auto native = Residual(state, history, alpha, time, data);
        Vector error(Size());
        for (int i = 0; i < system.rows.size(); ++i)
            error[system.rows[i]] = system.residual[i] - native[system.rows[i]];
        const auto constraints = Constraints(state, time, data);
        const auto norms       = Norms(error, constraints, m_robustUstar);
        emit("fixed_point", {{"block", J::Quote(name)},
                             {"scaled_difference", J::Numbers(norms)}});
        Require(
            Maximum(norms) <=
                1e-10 * std::max(1., Maximum(Norms(native, constraints,
                                                   m_robustUstar))),
            "frozen Picard block differs from the native same-state residual");
    };
    auto apply = [](Vector &state, const PicardBlock &block, const Vector &step,
                    double relaxation) {
        for (int i = 0; i < block.rows.size(); ++i)
            state[block.rows[i]] += relaxation * step[i];
    };
    auto scalarStep = [&](int b, Vector &state, double requested,
                          ScalarDomainResult &domain, std::string &ground) {
        const char *name = b == K() ? "k" : "epsilon_transport_trace";
        const auto block =
            AssemblePicardBlock(b, state, history, alpha, time, data, external);
        if (m_picardVerify && iterations == 0)
            blockCheck(block, state, name);
        auto rhs = block.residual;
        for (double &entry : rhs)
            entry = -entry;
        const auto step = Factor(block.matrix).Apply(rhs);
        auto defect     = Multiply(block.matrix, step);
        for (int i = 0; i < defect.size(); ++i)
            defect[i] -= rhs[i];
        Require(Maximum(defect) <= 1e-9 * std::max(1., Maximum(rhs)),
                "inaccurate Picard scalar block solve");
        ++(b == K() ? kSolves : epsilonSolves);
        emit("solve", {{"block", J::Quote(name)},
                       {"dimension", J::Number(step.size())},
                       {"frozen_scalar_defect", J::Number(Maximum(defect))}});
        for (int i = 0; i < 14; ++i, requested *= .5)
        {
            auto trial = state;
            apply(trial, block, step, requested);
            if (admissible(trial, domain, ground))
            {
                state = std::move(trial);
                return requested;
            }
        }
        emit("rejected_scalar",
             {{"block", J::Quote(name)}, {"domain", ScalarDomainJSON(domain)}});
        return 0.;
    };
    auto report = [&](const char *event, const Vector &norms,
                      const ScalarDomainResult &domain,
                      const std::string &ground, const Vector &relaxation) {
        double traceDefect = 0.;
        for (int id : m_traceIds)
        {
            const int row = offset[Epsilon()] + id;
            traceDefect   = std::max(traceDefect,
                                     std::abs(s[row] - domain.traceState[row]));
        }
        emit(event, {{"norms", J::Numbers(norms)},
                     {"relaxation", J::Numbers(relaxation)},
                     {"scalar_domain", ScalarDomainJSON(domain)},
                     {"ground", ground},
                     {"trace_coefficient_defect", J::Number(traceDefect)},
                     {"solve_counts",
                      J::Numbers(std::array<int, 4>{flowSolves, flowSolves,
                                                    kSolves, epsilonSolves})}});
    };
    auto perturb = [&](Vector &guess, double amplitude) {
        const std::array<double, 4> scales{
            {m_robustUstar, m_robustUstar,
             m_parameters.model.Q() * m_robustUstar * m_robustUstar,
             std::pow(m_robustUstar, 3) /
                 (m_parameters.model.Kappa() * m_parameters.geometry.height)}};
        for (int b = 0; b < InnerPressure(); ++b)
            for (int i = 0; i < size[b]; ++i)
            {
                const int row = offset[b] + i;
                if (!external.count(row))
                    guess[row] += amplitude * scales[b] * std::sin(1. + row);
            }
    };

    ScalarDomainResult domain;
    std::string ground;
    if (m_picardGuessPerturbation)
        perturb(s, m_picardGuessPerturbation);
    Require(admissible(s, domain, ground),
            "Picard initial state fails scalar or ground domain: " +
                domain.failure);
    m_lastFeasibleState = s;
    auto norms          = normsAt(s);
    auto historyValues  = history;
    for (auto &value : historyValues)
        if (value == 0.)
            value = 0.;
    emit("stage_input",
         {{"unknowns", J::Number(Size())},
          {"tolerance", J::Number(tolerance)},
          {"history_sha256", J::Quote(RANSWallR3::Hash(history))},
          {"history_value_sha256", J::Quote(RANSWallR3::Hash(historyValues))},
          {"incoming_guess_sha256", J::Quote(incomingGuessHash)},
          {"guess_sha256", J::Quote(RANSWallR3::Hash(s))}});
    report("initial", norms, domain, ground, {});
    if (m_picardVerify)
    {
        auto perturbed = s;
        perturb(perturbed, 1e-7);
        ScalarDomainResult perturbedDomain;
        std::string perturbedGround;
        Require(admissible(perturbed, perturbedDomain, perturbedGround),
                "Picard consistency perturbation leaves the scalar or ground "
                "domain");
        for (const int b : {0, K(), Epsilon()})
        {
            const auto baseName = b == 0     ? "initial_flow"
                                  : b == K() ? "initial_k"
                                             : "initial_epsilon";
            blockCheck(
                AssemblePicardBlock(b, s, history, alpha, time, data, external),
                s, baseName);
            const auto name = b == 0     ? "perturbed_flow"
                              : b == K() ? "perturbed_k"
                                         : "perturbed_epsilon";
            blockCheck(AssemblePicardBlock(b, perturbed, history, alpha, time,
                                           data, external),
                       perturbed, name);
        }
    }
    while (Maximum(norms) > tolerance && iterations < maximumIterations)
    {
        const auto flow =
            AssemblePicardBlock(0, s, history, alpha, time, data, external);
        const int nv = std::lower_bound(flow.rows.begin(), flow.rows.end(),
                                        DifferentialSize()) -
                       flow.rows.begin();
        const int np = flow.rows.size() - nv;
        WallMatrix a(nv, Vector(nv)), d(np, Vector(nv)), g(nv, Vector(np));
        for (int i = 0; i < nv; ++i)
        {
            std::copy_n(flow.matrix[i].begin(), nv, a[i].begin());
            std::copy_n(flow.matrix[i].begin() + nv, np, g[i].begin());
        }
        for (int i = 0; i < np; ++i)
            std::copy_n(flow.matrix[nv + i].begin(), nv, d[i].begin());
        Vector rhs(nv);
        for (int i = 0; i < nv; ++i)
            rhs[i] = -flow.residual[i];
        const Factor velocity(a);
        const auto predictor = velocity.Apply(rhs);
        emit("solve", {{"block", J::Quote("velocity_predictor")},
                       {"dimension", J::Number(nv)}});

        // A du + G dp = -r_u, D du + C dp = -r_p.
        // Use the same A inverse for both the Schur columns and correction.
        WallMatrix schur(np, Vector(np)), inverseGradient(np, Vector(nv));
        for (int j = 0; j < np; ++j)
        {
            for (int i = 0; i < nv; ++i)
                rhs[i] = g[i][j];
            inverseGradient[j] = velocity.Apply(rhs);
            const auto column  = Multiply(d, inverseGradient[j]);
            for (int i = 0; i < np; ++i)
                schur[i][j] = column[i] - flow.matrix[nv + i][nv + j];
        }
        auto pressureRhs = Multiply(d, predictor);
        for (int i = 0; i < np; ++i)
            pressureRhs[i] += flow.residual[nv + i];
        const auto pressure = Factor(schur).Apply(pressureRhs);
        Vector step(predictor);
        for (int i = 0; i < nv; ++i)
            for (int j = 0; j < np; ++j)
                step[i] -= inverseGradient[j][i] * pressure[j];
        step.insert(step.end(), pressure.begin(), pressure.end());
        ++flowSolves;
        auto frozenDefect = Multiply(flow.matrix, step);
        for (int i = 0; i < step.size(); ++i)
            frozenDefect[i] += flow.residual[i];
        Require(Maximum(frozenDefect) <=
                    1e-9 * std::max(1., Maximum(flow.residual)),
                "inaccurate Picard predictor/correction solve");
        emit("solve",
             {{"block", J::Quote("pressure_correction")},
              {"dimension", J::Number(np)},
              {"schur_velocity_rhs", J::Number(np)},
              {"predictor_continuity_gauge_max",
               J::Number(Maximum(pressureRhs))},
              {"frozen_flow_defect", J::Number(Maximum(frozenDefect))}});
        if (m_picardVerify && iterations == 0)
        {
            auto directRhs = flow.residual;
            for (double &entry : directRhs)
                entry = -entry;
            const auto direct = RANSWallAlgebra::Solve(flow.matrix, directRhs);
            Vector difference(step.size());
            double velocityCorrection = 0.;
            for (int i = 0; i < step.size(); ++i)
                difference[i] = step[i] - direct[i];
            for (int i = 0; i < nv; ++i)
                velocityCorrection = std::max(velocityCorrection,
                                              std::abs(step[i] - predictor[i]));
            auto predicted = s, corrected = s;
            apply(predicted, flow, predictor, 1.);
            apply(corrected, flow, step, 1.);
            auto pressureDirection = Vector(Size());
            for (int j = 0; j < np; ++j)
                pressureDirection[flow.rows[nv + j]] = pressure[j];
            const auto nativePressure =
                AffinePressureAction(pressureDirection, external);
            const auto assembledPressure = Multiply(g, pressure);
            Vector pressureError(nv);
            for (int i = 0; i < nv; ++i)
                pressureError[i] =
                    assembledPressure[i] - nativePressure[flow.rows[i]];
            const auto predictedNorms = normsAt(predicted),
                       correctedNorms = normsAt(corrected);
            auto directState          = s;
            apply(directState, flow, direct, 1.);
            const auto correctedInterface =
                Audit(corrected, history, alpha, time, data, m_robustUstar)
                    .at("interface");
            const auto directInterface =
                Audit(directState, history, alpha, time, data, m_robustUstar)
                    .at("interface");
            Vector interfaceError(correctedInterface.size());
            for (int i = 0; i < interfaceError.size(); ++i)
                interfaceError[i] = correctedInterface[i] - directInterface[i];
            emit("frozen_saddle_comparison",
                 {{"coefficient_difference", J::Number(Maximum(difference))},
                  {"velocity_correction", J::Number(velocityCorrection)},
                  {"affine_pressure_difference",
                   J::Number(Maximum(pressureError))},
                  {"interface_corrected", J::Numbers(correctedInterface)},
                  {"interface_direct", J::Numbers(directInterface)},
                  {"original_predictor_continuity_gauge",
                   J::Numbers(Vector(predictedNorms.begin() + InnerPressure(),
                                     predictedNorms.begin() + Gauge() + 1))},
                  {"original_corrected_continuity_gauge",
                   J::Numbers(Vector(correctedNorms.begin() + InnerPressure(),
                                     correctedNorms.begin() + Gauge() + 1))}});
            Require(Maximum(difference) <= 1e-8 * std::max(1., Maximum(direct)),
                    "Picard pressure correction disagrees with the direct "
                    "frozen saddle system");
            Require(Maximum(pressureError) <=
                        1e-10 * std::max(1., Maximum(assembledPressure)),
                    "Picard affine pressure increment disagrees with the "
                    "native action");
            Require(Maximum(interfaceError) <=
                        1e-8 * std::max(1., Maximum(directInterface)),
                    "Picard interface observations disagree with the direct "
                    "frozen saddle state");
        }

        const double merit =
            std::inner_product(norms.begin(), norms.end(), norms.begin(), 0.);
        bool accepted   = false;
        double fraction = 1.;
        for (int retry = 0; retry < 12; ++retry, fraction *= .5)
        {
            auto trial = s;
            ScalarDomainResult trialDomain;
            std::string trialGround;
            const double flowRelax = m_picardRelaxation[0] * fraction;
            apply(trial, flow, step, flowRelax);
            if (!admissible(trial, trialDomain, trialGround))
            {
                emit("rejected_flow",
                     {{"fraction", J::Number(fraction)},
                      {"domain", ScalarDomainJSON(trialDomain)}});
                continue;
            }
            const auto flowNorms = trialNorms(trial, fraction);
            if (flowNorms.empty())
                continue;
            const double kRelax =
                scalarStep(K(), trial, m_picardRelaxation[1] * fraction,
                           trialDomain, trialGround);
            if (!kRelax)
                continue;
            // The k update may leave a trace mismatch. Assemble epsilon
            // with the new k and freshly evaluated mechanical wall forcing
            // and sensors. No other unknown is lifted in this scalar solve.
            const double epsRelax =
                scalarStep(Epsilon(), trial, m_picardRelaxation[2] * fraction,
                           trialDomain, trialGround);
            if (!epsRelax)
                continue;
            const auto next = trialNorms(trial, fraction);
            if (next.empty())
                continue;
            const Vector relaxation{flowRelax, kRelax, epsRelax};
            emit("trial", {{"norms", J::Numbers(next)},
                           {"relaxation", J::Numbers(relaxation)},
                           {"post_flow_continuity_gauge",
                            J::Numbers(Vector{flowNorms[InnerPressure()],
                                              flowNorms[OuterPressure()],
                                              flowNorms[Gauge()]})}});
            if (Maximum(next) > tolerance &&
                std::inner_product(next.begin(), next.end(), next.begin(),
                                   0.) >= merit)
                continue;
            s                   = std::move(trial);
            m_lastFeasibleState = s;
            norms               = next;
            ++iterations;
            report("accepted_sweep", norms, trialDomain, trialGround,
                   relaxation);
            accepted = true;
            break;
        }
        Require(
            accepted,
            "Picard could not find an admissible residual-decreasing sweep");
    }
    Require(Maximum(norms) <= tolerance,
            "Picard sweep limit before full stage convergence");
    emit("converged", {{"norms", J::Numbers(norms)}});
    return s;
}
} // namespace Nektar
