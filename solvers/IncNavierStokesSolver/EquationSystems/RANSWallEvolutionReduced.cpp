///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallReducedModel.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <set>

namespace Nektar
{
using RANSWallAlgebra::Maximum;
using RANSWallAlgebra::Multiply;
using RANSWallAlgebra::Require;
namespace J = RANSWallSolveJson;
namespace
{
using Vector = RANSWallAlgebra::Vector;
double Dot(const Vector &a, const Vector &b)
{
    if (a.size() != b.size())
        throw std::logic_error("reduced dot product shape mismatch");
    long double sum = 0;
    for (int i = 0; i < a.size(); ++i)
        sum += (long double)a[i] * b[i];
    return static_cast<double>(sum);
}
double Length(const Vector &a)
{
    return std::sqrt(std::max(0., Dot(a, a)));
}
Vector HNorms(Vector n, int k)
{
    n.at(k) = 0.;
    return n;
}
void Apply(Vector &state, const std::vector<int> &rows, const Vector &step,
           double fraction = 1.)
{
    if (rows.size() != step.size())
        throw std::logic_error("reduced increment shape mismatch");
    for (int i = 0; i < rows.size(); ++i)
        state.at(rows[i]) += fraction * step[i];
}
} // namespace

RANSWallEvolutionPatch::PicardFlowResult RANSWallEvolutionPatch::
    PicardFlowCorrection(const Vector &state, const Vector &history,
                         double alpha, double time,
                         const RANSWallEvolutionData &data,
                         const std::map<int, double> &external) const
{
    PicardFlowResult result;
    result.block =
        AssemblePicardBlock(0, state, history, alpha, time, data, external);
    const auto &flow = result.block;
    const int nv     = std::lower_bound(flow.rows.begin(), flow.rows.end(),
                                        DifferentialSize()) -
                   flow.rows.begin();
    const int np        = flow.rows.size() - nv;
    result.velocitySize = nv;
    result.pressureSize = np;
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
    result.predictor = velocity.Apply(rhs);
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
    auto pressureRhs = Multiply(d, result.predictor);
    for (int i = 0; i < np; ++i)
        pressureRhs[i] += flow.residual[nv + i];
    result.predictorDefect = Maximum(pressureRhs);
    const auto pressure    = Factor(schur).Apply(pressureRhs);
    result.step            = result.predictor;
    for (int i = 0; i < nv; ++i)
        for (int j = 0; j < np; ++j)
            result.step[i] -= inverseGradient[j][i] * pressure[j];
    result.step.insert(result.step.end(), pressure.begin(), pressure.end());
    result.defect = Multiply(flow.matrix, result.step);
    for (int i = 0; i < result.step.size(); ++i)
        result.defect[i] += flow.residual[i];
    Require(Maximum(result.defect) <=
                1e-9 * std::max(1., Maximum(flow.residual)),
            "inaccurate reduced velocity predictor/correction solve");
    return result;
}

RANSWallEvolutionPatch::NonKResponse RANSWallEvolutionPatch::EquilibrateNonK(
    const Vector &state, const Vector &history, double alpha, double time,
    const RANSWallEvolutionData &data, double tolerance,
    int maximumUpdates) const
{
    Require(state.size() == Size() && history.size() == state.size() &&
                tolerance > 0. && std::isfinite(tolerance) &&
                maximumUpdates >= 0,
            "invalid controlled non-k response request");
    NonKResponse result;
    result.state = state;
    auto current = ExternalLift(state, time, data);
    Require(current == state, "controlled non-k input violates external data");
    const Vector held(current.begin() + offset[K()],
                      current.begin() + offset[Epsilon()]);
    auto external = ExternalConstraints(time, data);
    for (int id : m_traceIds)
        external.erase(offset[Epsilon()] + id);
    auto admissible = [&](const Vector &s) {
        bool ground = false;
        WallTraceJSON(s, ground);
        return ground && ScalarDomain(s, time, data).admissible;
    };
    auto evaluate = [&](const Vector &s) {
        result.residual = Residual(s, history, alpha, time, data);
        result.norms =
            Norms(result.residual, Constraints(s, time, data), m_robustUstar);
        result.error = Maximum(HNorms(result.norms, K()));
    };
    try
    {
        if (!admissible(current))
        {
            result.failure = "evaluation_or_ground_domain";
            return result;
        }
        evaluate(current);
        while (result.error > tolerance && result.updates < maximumUpdates)
        {
            const auto correction = PicardFlowCorrection(
                current, history, alpha, time, data, external);
            ++result.flowSolves;
            const double before =
                Dot(HNorms(result.norms, K()), HNorms(result.norms, K()));
            bool accepted = false;
            std::string trialFailure;
            for (int retry = 0; retry < 9; ++retry)
            {
                const double fraction = std::ldexp(1., -retry);
                auto trial            = current;
                Apply(trial, correction.block.rows, correction.step, fraction);
                try
                {
                    if (!admissible(trial))
                        continue;
                    const auto epsilon = AssemblePicardBlock(
                        Epsilon(), trial, history, alpha, time, data, external);
                    auto rhs = epsilon.residual;
                    for (auto &v : rhs)
                        v = -v;
                    const auto step = Factor(epsilon.matrix).Apply(rhs);
                    ++result.epsilonSolves;
                    Apply(trial, epsilon.rows, step, fraction);
                    if (!admissible(trial))
                        continue;
                    const auto residual =
                        Residual(trial, history, alpha, time, data);
                    const auto norms =
                        Norms(residual, Constraints(trial, time, data),
                              m_robustUstar);
                    const double error = Maximum(HNorms(norms, K()));
                    if (error > tolerance &&
                        Dot(HNorms(norms, K()), HNorms(norms, K())) >= before)
                        continue;
                    if (!std::equal(held.begin(), held.end(),
                                    trial.begin() + offset[K()]))
                        throw std::logic_error(
                            "controlled non-k solve modified k");
                    current         = std::move(trial);
                    result.residual = residual;
                    result.norms    = norms;
                    result.error    = error;
                    ++result.updates;
                    accepted = true;
                    break;
                }
                catch (const std::runtime_error &e)
                {
                    trialFailure = e.what();
                }
                catch (const std::invalid_argument &e)
                {
                    trialFailure = e.what();
                }
            }
            if (!accepted)
            {
                result.failure = "non_k_line_search" +
                                 (trialFailure.empty() ? std::string{}
                                                       : ": " + trialFailure);
                return result;
            }
        }
        if (result.error > tolerance)
        {
            result.failure = "non_k_iteration_limit";
            return result;
        }
        result.converged = true;
        result.state     = std::move(current);
    }
    catch (const std::runtime_error &e)
    {
        result.failure = std::string("non_k_numerical_failure: ") + e.what();
    }
    catch (const std::invalid_argument &e)
    {
        result.failure = std::string("non_k_domain_failure: ") + e.what();
    }
    return result;
}

RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::SolvePicardContinuation(
    Vector state, const Vector &history, double alpha, double time,
    const RANSWallEvolutionData &data, int &iterations, double tolerance,
    int maximumIterations) const
{
    const double eta = m_picardContinuationEta, operatorAlpha = alpha + eta;
    const double auxiliaryTolerance = eta > 0. && m_picardContinuationTolerance > 0.
                                          ? m_picardContinuationTolerance : tolerance;
    Require(std::isfinite(operatorAlpha), "nonfinite continued stage alpha");
    auto emit = [&](const std::string &event,
                    std::vector<std::pair<std::string,std::string>> fields) {
        fields.emplace_back("event", J::Quote(event));
        fields.emplace_back("physical_alpha", J::Number(alpha));
        fields.emplace_back("operator_alpha", J::Number(operatorAlpha));
        fields.emplace_back("eta", J::Number(eta));
        fields.emplace_back("original_tolerance", J::Number(tolerance));
        fields.emplace_back("auxiliary_tolerance", J::Number(auxiliaryTolerance));
        fields.emplace_back("time", J::Number(time));
        fields.emplace_back("history_sha256", J::Quote(RANSWallR3::Hash(history)));
        std::cout << "EP-AWF Continuation " << J::Object(fields) << std::endl;
    };
    emit("target", {{"physical_history_advanced", "false"}});
    auto root = SolveSegregatedPicard(std::move(state), history, operatorAlpha,
                                     time, data, iterations, auxiliaryTolerance, maximumIterations);
    m_lastFeasibleState = root;
    const auto constraints = Constraints(root, time, data);
    const auto original = Residual(root, history, alpha, time, data);
    const auto continued = Residual(root, history, operatorAlpha, time, data);
    const auto originalNorms = Norms(original, constraints, m_robustUstar);
    const auto continuedNorms = Norms(continued, constraints, m_robustUstar);
    Require(Maximum(continuedNorms) <= auxiliaryTolerance,
            "continued root missed its independently recomputed complete residual");
    Vector massRate(Size()), identity(Size());
    const auto lifted = ExternalLift(root, time, data);
    for (int b = 0; b < InnerPressure(); ++b)
    {
        Vector delta(size[b]);
        for (int i = 0; i < size[b]; ++i)
            delta[i] = lifted[offset[b]+i] - history[offset[b]+i];
        const auto load = Multiply(m_mass[b], delta);
        for (int i = 0; i < size[b]; ++i) massRate[offset[b]+i] = eta * load[i];
    }
    for (const auto &[row, value] : constraints) massRate[row] = 0.;
    for (int i = 0; i < Size(); ++i)
        identity[i] = continued[i] - original[i] - massRate[i];
    const double defect = Maximum(identity) /
        std::max({1., Maximum(original), Maximum(continued), Maximum(massRate)});
    Require(defect < 2e-12, "continued operator differs from the declared physical mass homotopy");
    const auto audit = CompleteRestrictedAudit(root, history, alpha, time, data, true);
    bool ground = false;
    const auto wall = WallTraceJSON(root, ground);
    Require(ground, "continued root is outside the continuous ground domain");
    std::ofstream out("picard-continuation-root.json");
    out << J::Object({{"status", J::Quote("algebraic_continuation_root_only")},
        {"state", J::Numbers(root)}, {"history", J::Numbers(history)},
        {"alpha", J::Number(alpha)}, {"operator_alpha", J::Number(operatorAlpha)},
        {"eta", J::Number(eta)}, {"time", J::Number(time)},
        {"require_history_h0", "false"}, {"inner_stress_model", J::Quote(InnerStressName())}}) << '\n';
    out.close();
    Require(out.good(), "failed to save algebraic continuation root");
    emit("root", {{"original_norms", J::Numbers(originalNorms)},
        {"continued_norms", J::Numbers(continuedNorms)},
        {"iterations", J::Number(iterations)}, {"mass_identity_defect", J::Number(defect)},
        {"ground", wall}, {"outer_range", J::Numbers(audit.at("r2_outer_range"))},
        {"gradient_witness", J::Numbers(audit.at("r4_outer_gradient_witness"))},
        {"qualified", J::Boolean(audit.at("r2_qualified")[0] == 1.)},
        {"scalar_domain", J::Numbers(audit.at("r2_scalar_domain"))},
        {"covariance", J::Numbers(audit.at("r2_covariance"))}});
    // Even an eta-zero qualified root from this diagnostic must be submitted
    // to the ordinary native stage path before physical history can advance.
    throw std::runtime_error("Picard algebraic continuation diagnostic completed a root without advancing physical history");
}

RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::SolveReducedPicard(
    Vector state, const Vector &history, double alpha, double time,
    const RANSWallEvolutionData &data, int &iterations, double tolerance,
    int maximumIterations) const
{
    const auto started = std::chrono::steady_clock::now();
    const std::string stageTag = "t" + J::Number(time) + "-" +
                                  RANSWallR3::Hash(history).substr(0,12);
    auto emit          = [&](const std::string &event,
                    std::vector<std::pair<std::string, std::string>> fields) {
        fields.emplace_back("event", J::Quote(event));
        fields.emplace_back("iteration", J::Number(iterations));
        fields.emplace_back("alpha", J::Number(alpha));
        fields.emplace_back("time", J::Number(time));
        fields.emplace_back(
            "seconds", J::Number(std::chrono::duration<double>(
                                     std::chrono::steady_clock::now() - started)
                                              .count()));
        std::cout << "EP-AWF Reduced " << J::Object(fields) << std::endl;
    };
    auto external = ExternalConstraints(time, data);
    for (int id : m_traceIds)
        external.erase(offset[Epsilon()] + id);
    std::vector<int> rows;
    for (int i = offset[K()]; i < offset[Epsilon()]; ++i)
        if (!external.count(i))
            rows.push_back(i);
    std::vector<int> index(size[K()], -1);
    for (int i = 0; i < rows.size(); ++i)
        index[rows[i] - offset[K()]] = i;
    const double kScale =
        m_parameters.model.Q() * m_robustUstar * m_robustUstar;
    const double residualScale = kScale * m_robustUstar /
                                 m_parameters.geometry.height *
                                 std::sqrt(m_volumes[K()]);
    const Factor mass(PicardKMass(rows));
    auto dualNorm = [&](const Vector &load) {
        return std::sqrt(std::max(0., Dot(load, mass.Apply(load)))) /
               residualScale;
    };
    auto extract = [&](const Vector &full) {
        Vector r;
        for (int row : rows)
            r.push_back(full.at(row));
        return r;
    };
    auto admissible = [&](const Vector &s) {
        bool ground = false;
        WallTraceJSON(s, ground);
        return ground && ScalarDomain(s, time, data).admissible;
    };
    int responseSolves = 0, flowSolves = 0, epsilonSolves = 0;
    struct Evaluation
    {
        bool valid = false;
        NonKResponse response;
        ScalarDomainResult domain;
        std::vector<OuterBound> bounds;
        Vector r;
        double merit = 1e300, violation = 1e300;
        std::array<double, 4> family{};
    };
    auto evaluate = [&](const Vector &guess, double accuracy, int budget,
                        const std::string &purpose) {
        Evaluation e;
        e.response = EquilibrateNonK(guess, history, alpha, time, data,
                                     accuracy, budget);
        ++responseSolves;
        flowSolves += e.response.flowSolves;
        epsilonSolves += e.response.epsilonSolves;
        emit("non_k_response",
             {{"purpose", J::Quote(purpose)},
              {"converged", J::Boolean(e.response.converged)},
              {"required", J::Number(accuracy)},
              {"attained", J::Number(e.response.error)},
              {"updates", J::Number(e.response.updates)},
              {"flow_solves", J::Number(e.response.flowSolves)},
              {"epsilon_solves", J::Number(e.response.epsilonSolves)},
              {"failure", J::Quote(e.response.failure)}});
        if (!e.response.converged)
        {
            if (e.response.state != guess)
                throw std::logic_error(
                    "failed non-k response changed caller state");
            return e;
        }
        e.domain = ScalarDomain(e.response.state, time, data, true);
        if (!e.domain.admissible)
            return e;
        bool ground = false;
        WallTraceJSON(e.response.state, ground);
        if (!ground)
            return e;
        e.bounds             = OuterBoundsR5(e.domain);
        const auto violation = OuterViolationR5(e.bounds);
        e.family             = violation.family;
        e.violation          = violation.Max();
        e.r                  = extract(e.response.residual);
        e.merit              = dualNorm(e.r);
        e.valid              = true;
        return e;
    };
    // Honor an explicitly tighter request, but do not automatically demand
    // H accuracy beyond the complete original stage tolerance at the end.
    const double minimumInnerAccuracy = std::min(tolerance, m_reducedHTolerance);
    double innerTolerance = m_reducedHTolerance;
    emit("protocol", {{"physical_history_sha256", J::Quote(RANSWallR3::Hash(history))},
        {"radius_scaling", J::Quote("Euclidean norm of free k coefficient displacement divided by reference K")},
        {"response_method", J::Quote("one-sided differences of converged sequential non-k solves")},
        {"half_response_warm_start", J::Quote("interpolate converged non-k coefficients; preserve exact requested k")},
        {"initial_basis_budget", J::Number(m_reducedBasis)},
        {"initial_restoration_weight", J::Number(m_reducedRestorationWeight)},
        {"response_krylov_budget", J::Number(m_reducedKrylov)},
        {"maximum_basis", J::Number(32)}, {"response_accuracy", J::Number(innerTolerance)}});
    auto base = evaluate(state, innerTolerance, m_reducedHUpdates, "base");
    Require(base.valid,
            "reduced base non-k solve failed: " + base.response.failure);
    state               = base.response.state;
    m_lastFeasibleState = state;
    auto report         = [&](const Evaluation &e, const std::string &event) {
        bool ground      = false;
        const auto wall  = WallTraceJSON(e.response.state, ground);
        const auto audit = CompleteRestrictedAudit(e.response.state, history,
                                                           alpha, time, data, true);
        if (m_reducedVerify)
        {
            const std::string file = "reduced-trial-" + stageTag + "-" + std::to_string(iterations) +
                                     "-" + event + ".json";
            std::ofstream out(file);
            out << J::Object({{"status", J::Quote("nonlinear_trial_only")},
                {"state", J::Numbers(e.response.state)}, {"history", J::Numbers(history)},
                {"alpha", J::Number(alpha)}, {"time", J::Number(time)},
                {"require_history_h0", "false"},
                {"inner_stress_model", J::Quote("SaturatedR2")}}) << '\n';
            out.close();
            Require(out.good(), "failed to save reduced nonlinear trial packet");
            emit("trial_packet", {{"path", J::Quote(file)}});
        }
        emit(event,
                     {{"norms", J::Numbers(e.response.norms)},
              {"reduced_norm", J::Number(e.merit)},
              {"violation", J::Number(e.violation)},
              {"families", J::Numbers(e.family)},
              {"qualification_inventory", OuterViolationR5(e.bounds).report},
              {"outer_range", J::Numbers(audit.at("r2_outer_range"))},
              {"covariance", J::Numbers(audit.at("r2_covariance"))},
              {"scalar_domain", ScalarDomainJSON(e.domain)},
              {"gradient_witness",
               J::Numbers(audit.at("r4_outer_gradient_witness"))},
              {"qualified", J::Boolean(audit.at("r2_qualified")[0] == 1.)},
              {"ground", wall},
              {"state_sha256", J::Quote(RANSWallR3::Hash(e.response.state))},
              {"work", J::Numbers(std::array<int, 3>{responseSolves, flowSolves,
                                                     epsilonSolves})}});
    };
    report(base, "initial");
    if (m_reducedVerify)
    {
        const auto repeated =
            evaluate(state, innerTolerance, m_reducedHUpdates, "repeat_same_k");
        Require(repeated.valid, "repeat reduced evaluation failed");
        auto discrepancy = repeated.r;
        for (int i = 0; i < rows.size(); ++i)
            discrepancy[i] -= base.r[i];
        auto disturbed = state;
        for (int i = 0; i < offset[K()]; ++i)
            if (!external.count(i))
                disturbed[i] += 1e-6 * m_robustUstar * std::sin(1. + i);
        const auto failed = EquilibrateNonK(disturbed, history, alpha, time,
                                            data, innerTolerance, 0);
        Require(!failed.converged && failed.state == disturbed &&
                    failed.failure == "non_k_iteration_limit",
                "insufficient response accuracy was silently accepted");
        const auto tightened =
            evaluate(disturbed, std::max(minimumInnerAccuracy, innerTolerance * .05),
                     m_reducedHUpdates, "tightened_same_k");
        Require(tightened.valid, "tight repeated reduced evaluation failed");
        auto tightDifference = tightened.r;
        for (int i = 0; i < rows.size(); ++i)
            tightDifference[i] -= base.r[i];
        emit("repeatability",
             {{"same_k_difference", J::Number(dualNorm(discrepancy))},
              {"tightened_k_difference", J::Number(dualNorm(tightDifference))},
              {"failed_response_rolled_back", "true"},
              {"k_coefficients_fixed", "true"}});
    }
    struct Column
    {
        std::string name;
        Vector k, r, observations, stateAction, traceAction;
        Evaluation shifted;
        double h = 0., error = 0.;
        double observationError = 0.;
    };
    struct Cut
    {
        Vector a;
        double b;
        int observation = -1, sign = 0;
        std::string kind;
    };
    double radius = m_reducedRadius;
    std::vector<std::pair<double, double>> filter{{base.merit, base.violation}};
    Vector previous;
    int failedModels = 0;
    int basisBudget = m_reducedBasis;
    double restorationWeight = m_reducedRestorationWeight;
    int qualificationStalls = 0;
    std::vector<OuterBound> activeWitnesses;
    auto qualifiedRoot = [&]() {
        if (Maximum(base.response.norms) > tolerance) return false;
        const auto audit = CompleteRestrictedAudit(state, history, alpha, time, data, true);
        return audit.at("r2_qualified")[0] == 1.;
    };
    while (iterations < maximumIterations)
    {
        // Original equations and qualification are the acceptance authority,
        // before any optional tightening of the conditional response accuracy.
        if (qualifiedRoot())
        {
            report(base, "converged");
            return state;
        }
        if (base.violation > 0.)
            innerTolerance = std::max(minimumInnerAccuracy,
                std::min(innerTolerance, 1e-7 * base.violation));
        if (base.response.error > innerTolerance)
        {
            auto refined = evaluate(state, innerTolerance, m_reducedHUpdates,
                                    "refine_base_accuracy");
            Require(refined.valid, "reduced base accuracy refinement failed: " +
                                       refined.response.failure);
            base = std::move(refined);
            state = base.response.state;
            m_lastFeasibleState = state;
            if (qualifiedRoot())
            {
                report(base, "converged");
                return state;
            }
        }
        const bool restore  = base.violation > 0.;
        const int dimension = rows.size();
        const auto contacts = WallTraceR4(state);
        std::vector<std::pair<std::string, Vector>> seeds;
        auto add = [&](const std::string &name, Vector v) {
            const double length = Length(v);
            if (length > 0.)
            {
                for (double &x : v)
                    x /= length;
                seeds.emplace_back(name, std::move(v));
            }
        };
        auto rhs = base.r;
        for (double &x : rhs)
            x = -x;
        WallMatrix conditionalJacobian;
        std::unique_ptr<Factor> responsePreconditioner;
        bool conditionalPreconditioner = false;
        try
        {
            conditionalJacobian = PicardKJacobian(state, alpha, time, data, rows);
            auto factor = std::make_unique<Factor>(conditionalJacobian);
            add("conditional_newton", factor->Apply(rhs));
            if (m_reducedKrylov)
            {
                responsePreconditioner = std::move(factor);
                conditionalPreconditioner = true;
            }
        }
        catch (const std::runtime_error &e)
        {
            emit("conditional_seed_failed", {{"failure", J::Quote(e.what())}});
        }
        const auto frozen = AssemblePicardBlock(K(), state, history, alpha,
                                                time, data, external);
        try
        {
            auto factor = std::make_unique<Factor>(frozen.matrix);
            add("rate_picard", factor->Apply(rhs));
            if (m_reducedKrylov && !responsePreconditioner)
                responsePreconditioner = std::move(factor);
        }
        catch (const std::runtime_error &e)
        {
            emit("frozen_seed_failed", {{"failure", J::Quote(e.what())}});
        }
        Vector lift(dimension);
        const double left = lower[0], right = upper[0], length = right - left;
        const auto groundData = fields[K()]->Boundary(
            [&](const RANSWallCoordinate &p) {
                return (p[0] - left) * (right - p[0]) / (length * length);
            },
            [&](int axis, int side, const RANSWallGlobalField::Element &e) {
                return axis == 2 && !side && e.lower[2] == 0.;
            });
        for (auto [id, value] : groundData)
            if (index[id] >= 0)
                lift[index[id]] = value;
        add("ground_lift", lift);
        auto pointSeed = [&](double x, double z, const std::string &name,
                             bool derivative = false) {
            const auto q = Point(K(), {x, 0., z},
                                 z == m_parameters.geometry.height ? 1 : 0);
            Vector v(dimension);
            for (int j = 0; j < q.ids.size(); ++j)
                if (index[q.ids[j]] >= 0)
                    v[index[q.ids[j]]] =
                        derivative ? q.derivative[0][j] : q.phi[j];
            add(name, std::move(v));
        };
        for (double x : {-2.6, 1.6})
        {
            if (x <= left || x >= right)
                continue;
            pointSeed(x, m_parameters.geometry.sampleHeight,
                      "opposite_sample_" + std::to_string(x));
            pointSeed(x, .04, "normalization_" + std::to_string(x));
        }
        // Include spatially distinct signed observations, not only the maximum.
        std::vector<int> ordered(base.bounds.size());
        std::iota(ordered.begin(), ordered.end(), 0);
        auto excess = [](const OuterBound &b) {
            return std::max((b.value - b.upper) / b.width,
                            (b.lower - b.value) / b.width);
        };
        std::stable_sort(ordered.begin(), ordered.end(), [&](int a, int b) {
            return excess(base.bounds[a]) > excess(base.bounds[b]);
        });
        const auto epsilonWitness = std::find_if(ordered.begin(), ordered.end(),
            [&](int id) { return base.bounds[id].family == 2; });
        if (epsilonWitness != ordered.end())
        {
            const double x = base.bounds[*epsilonWitness].row.coordinate[0];
            pointSeed(x, m_parameters.geometry.sampleHeight, "epsilon_lower_sample");
            pointSeed(x, m_parameters.geometry.height, "epsilon_lower_trace");
        }
        if (basisBudget >= 20 || failedModels)
        {
            auto dual = mass.Apply(base.r);
            Vector gradient(dimension);
            if (!conditionalJacobian.empty())
            {
                for (int i = 0; i < dimension; ++i)
                    for (int j = 0; j < dimension; ++j)
                        gradient[i] -= conditionalJacobian[j][i] * dual[j];
                add("conditional_mass_gradient", mass.Apply(gradient));
            }
            for (double &value : dual) value = -value;
            add("mass_residual", dual);
            std::vector<int> defects(dimension);
            std::iota(defects.begin(), defects.end(), 0);
            const int defectCount = std::min(2, dimension);
            std::partial_sort(defects.begin(), defects.begin()+defectCount, defects.end(),
                [&](int a, int b) { return std::abs(dual[a]) > std::abs(dual[b]); });
            for (int j = 0; j < defectCount; ++j)
            {
                const int i = defects[j];
                Vector mode(dimension); mode[i] = dual[i] < 0. ? -1. : 1.;
                add("modal_defect_" + std::to_string(rows[i]), std::move(mode));
            }
            auto ground = contacts;
            std::stable_sort(ground.begin(), ground.end(), [](const auto &a, const auto &b) {
                return a.value - a.error < b.value - b.error;
            });
            std::vector<double> locations;
            for (const auto &p : ground)
            {
                const double x = p.row.coordinate[0];
                if (x <= left || x >= right || p.value > .05*kScale ||
                    std::any_of(locations.begin(), locations.end(), [&](double y) { return std::abs(x-y)<.1; }))
                    continue;
                pointSeed(x, 0., "contact_value_" + std::to_string(locations.size()));
                pointSeed(x, 0., "contact_tangent_" + std::to_string(locations.size()), true);
                const double nearby = x + .02*length < right ? x + .02*length : x - .02*length;
                pointSeed(nearby, 0., "contact_nearby_" + std::to_string(locations.size()));
                locations.push_back(x);
                if (locations.size() == 2) break;
            }
        }
        int activeSeed = 0;
        for (const auto &b : activeWitnesses)
        {
            const double x = b.row.coordinate[0];
            pointSeed(x, m_parameters.geometry.sampleHeight,
                      "active_sample_" + std::to_string(activeSeed), b.family == 0);
            pointSeed(x, b.row.coordinate[2],
                      "active_location_" + std::to_string(activeSeed));
            ++activeSeed;
        }
        std::set<std::pair<int, int>> patches;
        int witness = 0;
        for (int id : ordered)
        {
            const auto &b = base.bounds[id];
            const int cell =
                static_cast<int>(std::floor((b.row.coordinate[0] - left) * 4));
            if (!patches.emplace(b.family, cell).second)
                continue;
            const double x = b.row.coordinate[0];
            pointSeed(x, m_parameters.geometry.sampleHeight,
                      "sample_witness_" + std::to_string(witness),
                      b.family == 0);
            pointSeed(x, m_parameters.geometry.height,
                      "trace_witness_" + std::to_string(witness));
            if (++witness == std::max(5, (basisBudget+1)/2))
                break;
        }
        if (!previous.empty()) add("previous", previous);
        std::vector<Vector> orthogonal;
        std::vector<Column> columns;
        double maxActionError = 0.;
        int generatedKrylov = 0;
        for (int seedId = 0; seedId < seeds.size(); ++seedId)
        {
            auto [name, seed] = seeds[seedId];
            // A positive physical ground polynomial supplies feasible one-sided
            // probes.
            double inward = 0.;
            for (const auto &contact : contacts)
                if (contact.value < 1e-6 * kScale)
                {
                    double a = 0., l = 0.;
                    for (int j = 0; j < contact.row.ids.size(); ++j)
                    {
                        const int i = index[contact.row.ids[j]];
                        if (i >= 0)
                        {
                            a += contact.row.phi[j] * seed[i];
                            l += contact.row.phi[j] * lift[i];
                        }
                    }
                    if (l > 1e-12)
                        inward = std::max(inward, (-a + .002) / l);
                }
            for (int i = 0; i < dimension; ++i)
                seed[i] += inward * lift[i];
            const double length = Length(seed);
            if (length < 1e-14)
                continue;
            for (double &x : seed)
                x /= length;
            auto q = seed;
            for (int pass = 0; pass < 2; ++pass)
                for (const auto &v : orthogonal)
                {
                    const double coefficient = Dot(q, v);
                    for (int i = 0; i < dimension; ++i)
                        q[i] -= coefficient * v[i];
                }
            if (Length(q) < 1e-7)
                continue;
            const double qLength = Length(q);
            for (double &x : q)
                x /= qLength;
            Column column;
            column.name = name;
            column.k    = seed;
            Evaluation shifted;
            double h = m_reducedStep;
            for (int retry = 0; retry < 5; ++retry, h *= .5)
            {
                auto guess = state;
                Apply(guess, rows, seed, kScale * h);
                if (!admissible(guess))
                    continue;
                shifted = evaluate(guess, innerTolerance, m_reducedHUpdates,
                                   "response_" + name);
                if (shifted.valid)
                    break;
            }
            if (!shifted.valid)
                continue;
            column.h       = h;
            column.shifted = std::move(shifted);
            column.r       = column.shifted.r;
            for (int i = 0; i < dimension; ++i)
                column.r[i] = (column.r[i] - base.r[i]) / h;
            Vector ds = column.shifted.response.state,
                   dt = column.shifted.domain.traceState;
            for (int i = 0; i < Size(); ++i)
            {
                ds[i] -= state[i];
                dt[i] -= base.domain.traceState[i];
            }
            column.observations.reserve(base.bounds.size());
            for (const auto &b : base.bounds)
                column.observations.push_back(
                    b.Action(b.projected ? dt : ds, offset[b.block]) / h);
            // Every direction used by the model receives a second, tighter
            // one-sided check; verification controls extra repeatability checks.
            {
                auto guess = state;
                Apply(guess, rows, seed, kScale * h * .5);
                const auto canonicalGuess = guess;
                // The full response is already on the local H branch. Its
                // midpoint predicts the non-k half response to O(h^2), while
                // the independently prescribed k coefficients remain exact.
                // This only changes the warm start: H must still converge to
                // the same declared tighter accuracy before the action is used.
                for (int j = 0; j < Size(); ++j)
                    if (j < offset[K()] || j >= offset[Epsilon()])
                        guess[j] += .5 * (column.shifted.response.state[j] - state[j]);
                const double halfAccuracy =
                    std::max(minimumInnerAccuracy, innerTolerance * .05);
                auto half = evaluate(
                    guess, halfAccuracy,
                    m_reducedHUpdates, "half_response_" + name);
                if (!half.valid)
                    half = evaluate(canonicalGuess, halfAccuracy,
                                    m_reducedHUpdates, "retry_half_response_" + name);
                if (!half.valid)
                    continue;
                // Independently check the accelerated warm start on the first
                // usable column of each verified model. Failure to recover the
                // same response branch rejects this column, rather than mixing
                // different local Y(k) solutions in a reduced model.
                if (m_reducedVerify && columns.empty())
                {
                    const auto canonical = evaluate(
                        canonicalGuess, halfAccuracy, m_reducedHUpdates,
                        "canonical_half_response_" + name);
                    if (!canonical.valid)
                        continue;
                    auto difference = half.r;
                    auto change = half.r;
                    for (int j = 0; j < dimension; ++j)
                    {
                        difference[j] = (difference[j] - canonical.r[j]) / (h * .5);
                        change[j] = (change[j] - base.r[j]) / (h * .5);
                    }
                    const double actionDifference = dualNorm(difference) /
                                                    std::max(1., dualNorm(change));
                    double observationDifference = 0.;
                    Vector ds = half.response.state, dt = half.domain.traceState;
                    for (int j = 0; j < Size(); ++j)
                    {
                        ds[j] -= canonical.response.state[j];
                        dt[j] -= canonical.domain.traceState[j];
                    }
                    for (const auto &b : base.bounds)
                        observationDifference = std::max(observationDifference,
                            std::abs(b.Action(b.projected ? dt : ds,
                                              offset[b.block])) / (h * .5));
                    const bool sameBranch = actionDifference <= .01 &&
                                            observationDifference <= .01;
                    emit("half_response_warm_start",
                        {{"name", J::Quote(name)},
                         {"required", J::Number(halfAccuracy)},
                         {"attained", J::Number(half.response.error)},
                         {"canonical_attained", J::Number(canonical.response.error)},
                         {"updates", J::Number(half.response.updates)},
                         {"canonical_updates", J::Number(canonical.response.updates)},
                         {"relative_action_difference", J::Number(actionDifference)},
                         {"observation_action_difference", J::Number(observationDifference)},
                         {"same_local_response", J::Boolean(sameBranch)}});
                    if (!sameBranch)
                        continue;
                }
                auto action = half.r;
                for (int i = 0; i < dimension; ++i)
                    action[i] = (action[i] - base.r[i]) / (h * .5);
                auto difference = action;
                for (int i = 0; i < dimension; ++i)
                    difference[i] -= column.r[i];
                column.error =
                    dualNorm(difference) / std::max(1., dualNorm(action));
                maxActionError = std::max(maxActionError, column.error);
                Vector hs = half.response.state, ht = half.domain.traceState;
                for (int i = 0; i < Size(); ++i)
                {
                    hs[i] -= state[i];
                    ht[i] -= base.domain.traceState[i];
                }
                double observationError = 0.;
                for (int j = 0; j < base.bounds.size(); ++j)
                {
                    const auto &b = base.bounds[j];
                    const double observed =
                        b.Action(b.projected ? ht : hs, offset[b.block]) /
                        (h * .5);
                    observationError =
                        std::max(observationError,
                                 std::abs(observed - column.observations[j]) /
                                     std::max(1., std::abs(observed)));
                }
                std::vector<std::string> observationChecks;
                for (int family = 0; family < 4; ++family)
                    for (int sign : {-1, 1})
                    {
                        int witness = -1;
                        double worst = -1e300;
                        for (int j = 0; j < base.bounds.size(); ++j)
                        {
                            const auto &b = base.bounds[j];
                            const double excess = sign * (b.value -
                                (sign < 0 ? b.lower : b.upper)) / b.width;
                            if (b.family == family && excess > worst)
                            {
                                witness = j;
                                worst = excess;
                            }
                        }
                        if (witness < 0) continue;
                        const auto &b = base.bounds[witness];
                        const double measured = b.Action(b.projected ? ht : hs,
                                                          offset[b.block]);
                        observationChecks.push_back(J::Object({
                            {"family", J::Number(family)}, {"sign", J::Number(sign)},
                            {"point", J::Numbers(b.row.coordinate)},
                            {"projected", J::Boolean(b.projected)},
                            {"signed_base", J::Number(worst)},
                            {"predicted_change", J::Number(sign * column.observations[witness] * h * .5)},
                            {"measured_change", J::Number(sign * measured)}}));
                    }
                Vector conditionalAction = conditionalJacobian.empty() ? Vector(dimension)
                    : Multiply(conditionalJacobian, seed);
                auto frozenAction = Multiply(frozen.matrix, seed);
                auto feedback = action;
                for (int i = 0; i < dimension; ++i)
                {
                    conditionalAction[i] *= kScale;
                    frozenAction[i] *= kScale;
                    feedback[i] -= conditionalAction[i];
                }
                emit("response_check",
                     {{"name", J::Quote(name)},
                      {"h", J::Number(h)},
                      {"relative_error", J::Number(column.error)},
                      {"conditional_action_norm", J::Number(dualNorm(conditionalAction))},
                      {"conditional_action_available", J::Boolean(!conditionalJacobian.empty())},
                      {"frozen_action_norm", J::Number(dualNorm(frozenAction))},
                      {"reduced_action_norm", J::Number(dualNorm(action))},
                      {"feedback_action_norm", J::Number(dualNorm(feedback))},
                      {"conditional_merit_slope", J::Number(Dot(base.r, mass.Apply(conditionalAction)) /
                                                             (residualScale * residualScale))},
                      {"signed_merit_slope",
                       J::Number(Dot(base.r, mass.Apply(column.r)) /
                                 (residualScale * residualScale))},
                      {"actual_merit_slope",
                       J::Number((column.shifted.merit * column.shifted.merit -
                                  base.merit * base.merit) / (2. * h))},
                      {"half_merit_slope",
                       J::Number((half.merit * half.merit - base.merit * base.merit) / h)},
                      {"observation_relative_error",
                       J::Number(observationError)},
                      {"usable", J::Boolean(column.error <= .03 && observationError <= .05)},
                      {"signed_observations", J::List(observationChecks)},
                      {"base_violation", J::Number(base.violation)},
                      {"shift_violation", J::Number(column.shifted.violation)},
                      {"half_violation", J::Number(half.violation)},
                      {"shift_domain", ScalarDomainJSON(column.shifted.domain)},
                      {"half_domain", ScalarDomainJSON(half.domain)},
                      {"profile_branch_change",
                       J::Boolean(column.shifted.domain.profileBranches !=
                                  base.domain.profileBranches)}});
                if (column.error > .03 || observationError > .05)
                    continue;
                column.observationError = observationError;
                column.r = std::move(action);
                for (int j = 0; j < base.bounds.size(); ++j)
                {
                    const auto &b = base.bounds[j];
                    column.observations[j] = b.Action(b.projected ? ht : hs,
                                                       offset[b.block]) / (h * .5);
                }
                column.h = h * .5;
                column.shifted = std::move(half);
            }
            column.stateAction = column.shifted.domain.effectiveState;
            column.traceAction = column.shifted.domain.traceState;
            for (int j = 0; j < Size(); ++j)
            {
                column.stateAction[j] = (column.stateAction[j] - base.domain.effectiveState[j]) / column.h;
                column.traceAction[j] = (column.traceAction[j] - base.domain.traceState[j]) / column.h;
            }
            // Retain the converged coefficients and actions, not a complete
            // duplicate of the large qualification inventory for every column.
            column.shifted.bounds.clear();
            column.shifted.bounds.shrink_to_fit();
            column.shifted.domain.observations.clear();
            column.shifted.domain.observations.shrink_to_fit();
            // Enrich beyond point/witness rows with the measured response
            // missing from the restricted k preconditioner. For a smooth
            // response this generates a preconditioned Krylov space for Jred.
            // It uses no all-field matrix; every new seed still receives its
            // own feasible, converged, two-scale response check above.
            if (responsePreconditioner && generatedKrylov < m_reducedKrylov &&
                (columns.empty() || name.rfind("response_krylov_", 0) == 0))
            {
                auto feedback = Multiply(conditionalPreconditioner ?
                    conditionalJacobian : frozen.matrix, column.k);
                for (int j = 0; j < dimension; ++j)
                    feedback[j] = column.r[j] - kScale * feedback[j];
                const double feedbackNorm = dualNorm(feedback);
                const double uncertainty = column.error * std::max(1., dualNorm(column.r));
                if (feedbackNorm > 10. * uncertainty)
                {
                    try
                    {
                        auto next = responsePreconditioner->Apply(feedback);
                        for (double &v : next) v /= -kScale;
                        const double norm = Length(next);
                        if (norm > 1e-14)
                        {
                            for (double &v : next) v /= norm;
                            const auto nextName = "response_krylov_" + std::to_string(++generatedKrylov);
                            seeds.insert(seeds.begin()+seedId+1, {nextName, std::move(next)});
                            emit("response_krylov_seed", {{"name", J::Quote(nextName)},
                                {"parent", J::Quote(name)},
                                {"preconditioner", J::Quote(conditionalPreconditioner ? "analytic_conditional_k" : "rate_picard_k")},
                                {"missing_response_norm", J::Number(feedbackNorm)},
                                {"two_scale_error_estimate", J::Number(uncertainty)}});
                        }
                    }
                    catch (const std::runtime_error &e)
                    {
                        emit("response_krylov_failed", {{"failure", J::Quote(e.what())}});
                    }
                }
            }
            orthogonal.push_back(std::move(q));
            columns.push_back(std::move(column));
            if (columns.size() >=
                std::min(32, basisBudget + failedModels * 3))
                break;
        }
        if (columns.size() < 3)
        {
            ++failedModels;
            innerTolerance = std::max(minimumInnerAccuracy, innerTolerance * .1);
            emit("inaccurate_or_failed_responses", {{"usable", J::Number(columns.size())},
                 {"required_accuracy", J::Number(innerTolerance)}});
            Require(failedModels < 4,
                    "insufficient accurate independent converged reduced responses");
            continue;
        }
        const int n = columns.size(), variables = n + (restore ? 1 : 0);
        // Do not spend successive response models approaching an exact
        // qualification boundary from outside at the scale of H-solve error.
        // This is an inward search target only. It vanishes with the original
        // residual, and neither the actual inventory nor final acceptance is
        // changed. Inside the envelope the original model bounds are used.
        const double qualificationBuffer = restore && base.violation < 1e-4
            ? std::min(1e-6, 1e-4 * base.merit) : 0.;
        // Keep both goals outside the envelope. Normalize the residual term
        // so its units do not overwhelm the signed qualification slack.
        const double residualWeight = restore ? restorationWeight * std::pow(.1, failedModels) /
                                                    std::max(1., base.merit * base.merit) : 1.;
        WallMatrix gram(n, Vector(n)), hessian(variables, Vector(variables));
        Vector gradient(variables), dualResidual = mass.Apply(base.r);
        std::vector<Vector> actionOrthogonal;
        for (int i = 0; i < n; ++i)
        {
            const auto dual = mass.Apply(columns[i].r);
            gradient[i]     = residualWeight * Dot(columns[i].r, dualResidual) /
                                        (residualScale * residualScale);
            auto q          = columns[i].r;
            for (int pass = 0; pass < 2; ++pass)
                for (const auto &v : actionOrthogonal)
                {
                    const double a = Dot(q, mass.Apply(v));
                    for (int j = 0; j < dimension; ++j)
                        q[j] -= a * v[j];
                }
            const double length =
                std::sqrt(std::max(0., Dot(q, mass.Apply(q))));
            if (length > 1e-7 * std::max(1., std::sqrt(std::max(
                                                 0., Dot(columns[i].r, dual)))))
            {
                for (double &x : q)
                    x /= length;
                actionOrthogonal.push_back(std::move(q));
            }
            for (int j = 0; j < n; ++j)
            {
                gram[i][j] = Dot(columns[i].k, columns[j].k);
                const double objective =
                    Dot(columns[j].r, dual) / (residualScale * residualScale);
                hessian[i][j] = residualWeight * objective + 1e-3 * gram[i][j]
                              + (i == j ? 1e-10 : 0.);
            }
        }
        if (restore)
            hessian[n][n] = 1.;
        std::vector<Cut> cuts;
        std::map<int, WallExtremum> groundLocations;
        std::map<int, ScalarObservation> domainLocations;
        auto insertCut = [&](Vector a, double b, int observation, int sign,
                             const std::string &kind) {
            const double norm = Length(a);
            if (norm < 1e-15)
                return;
            for (double &x : a)
                x /= norm;
            cuts.push_back({std::move(a), b / norm, observation, sign, kind});
        };
        for (int j = 0; j < base.bounds.size(); ++j)
        {
            const auto &b = base.bounds[j];
            for (int sign : {-1, 1})
            {
                Vector a(variables);
                for (int i = 0; i < n; ++i)
                    a[i] = sign * columns[i].observations[j];
                if (restore)
                    a[n] = -1.;
                insertCut(std::move(a),
                          -sign * (b.value - (sign < 0 ? b.lower : b.upper)) /
                              b.width - qualificationBuffer,
                          j, sign, "qualification");
            }
        }
        // Linear positivity cuts at the entire current scalar inventory.
        // For moving-ray observations these constrain the current locations;
        // every actual candidate also rebuilds and scans its moving inventory.
        // No positive floor is introduced for strict observations: strictness
        // is checked by the native evaluation-domain scan before evaluation.
        for (const auto &p : base.domain.observations)
        {
            const auto row = Point(p.block, p.coordinate, p.side);
            Vector a(variables);
            for (int i = 0; i < n; ++i)
                a[i] = -row.Value(p.dependentTrace ? columns[i].traceAction
                                                  : columns[i].stateAction,
                                  offset[p.block]) / p.scale;
            const int before = cuts.size();
            insertCut(std::move(a), p.value / p.scale, -1, 0, "evaluation_domain");
            if (cuts.size() > before) domainLocations.emplace(before, p);
        }
        if (restore)
        {
            Vector a(variables);
            a[n] = -1.;
            insertCut(std::move(a), 0., -1, 0, "restoration_slack");
            // Match the actual outside-envelope filter during construction,
            // rather than reject a useful residual direction only because
            // its full unconstrained proposal exceeds the allowed increase.
            Vector ceiling(variables);
            ceiling[n] = 1.;
            insertCut(std::move(ceiling), base.violation * 1.01 + 1e-7 + qualificationBuffer,
                      -1, 0, "qualification_funnel");
        }
        auto groundCut = [&](const WallExtremum &p) {
            Vector a(variables);
            for (int i = 0; i < n; ++i)
                for (int j = 0; j < p.row.ids.size(); ++j)
                {
                    const int c = index[p.row.ids[j]];
                    if (c >= 0)
                        a[i] -= kScale * p.row.phi[j] * columns[i].k[c];
                }
            const int before = cuts.size();
            insertCut(std::move(a), p.value - p.error, -1, 0,
                      "continuous_ground");
            if (cuts.size() > before) groundLocations.emplace(before,p);
        };
        for (const auto &p : contacts)
            groundCut(p);
        Vector solution(variables);
        if (restore)
            solution[n] = base.violation + qualificationBuffer + 1e-10;
        std::vector<int> active;
        Vector multipliers;
        bool qpOK   = false;
        bool groundModelOK = false;
        int qpSteps = 0, groundAdded = 0;
        RANSWallReduced::QuadraticResult quadratic;
        auto solveQP = [&]() {
            quadratic = RANSWallReduced::SolveQuadratic<Factor>(
                hessian, gradient, cuts, solution, 500);
            qpSteps += quadratic.steps;
            solution = quadratic.solution;
            active = quadratic.active;
            multipliers = quadratic.multipliers;
            return quadratic.solved;
        };
        WallMatrix kColumns;
        for (const auto &column : columns) kColumns.push_back(column.k);
        Vector displacement(dimension);
        for (int cutLoop = 0; cutLoop < 32; ++cutLoop)
        {
            try
            {
                qpOK = solveQP();
            }
            catch (const std::runtime_error &error)
            {
                emit("qp_numerical_failure", {{"failure", J::Quote(error.what())}});
                qpOK = false;
            }
            if (!qpOK)
                break;
            displacement = RANSWallReduced::BoundCorrection(solution, kColumns, radius);
            auto guess = state;
            Apply(guess, rows, displacement, kScale);
            const auto extrema = WallTraceR4(guess);
            const auto worst   = std::min_element(
                extrema.begin(), extrema.end(),
                [](const auto &a, const auto &b) {
                    return a.value - a.error < b.value - b.error;
                });
            if (worst->value >= worst->error)
            {
                groundModelOK = true;
                break;
            }
            const double roundoff = 256. * std::numeric_limits<double>::epsilon() *
                std::max(kScale, Maximum(extract(guess)));
            emit("ground_cut", {{"minimum", J::Number(worst->value)},
                 {"arithmetic_scale", J::Number(worst->error)},
                 {"point", J::Numbers(worst->row.coordinate)},
                 {"cut_loop", J::Number(cutLoop)},
                 {"double_roundoff_scale", J::Number(roundoff)}});
            if (worst->error - worst->value <= roundoff)
            {
                const auto liftColumn = std::find_if(columns.begin(), columns.end(),
                    [](const Column &c) { return c.name == "ground_lift"; });
                if (liftColumn != columns.end())
                {
                    const int c = liftColumn - columns.begin();
                    double inwardSlope = 0.;
                    for (int j = 0; j < worst->row.ids.size(); ++j)
                    {
                        const int i = index[worst->row.ids[j]];
                        if (i >= 0) inwardSlope += kScale * worst->row.phi[j] * columns[c].k[i];
                    }
                    if (inwardSlope > 0.)
                    {
                        double tau = std::max(2. * static_cast<double>(worst->error - worst->value) / inwardSlope,
                            8. * std::numeric_limits<double>::epsilon() * std::max(1., std::abs(solution[c])));
                        for (int retry = 0; retry < 12; ++retry, tau *= 2.)
                        {
                            auto coefficients = solution;
                            coefficients[c] += tau;
                            auto adjusted = RANSWallReduced::BoundCorrection(coefficients, kColumns, radius);
                            auto interior = state;
                            Apply(interior, rows, adjusted, kScale);
                            bool valid = false;
                            const auto wall = WallTraceJSON(interior, valid);
                            if (!valid) continue;
                            solution = std::move(coefficients);
                            displacement = std::move(adjusted);
                            groundModelOK = true;
                            emit("ground_inward_roundoff", {{"coefficient", J::Number(tau)},
                                 {"before_minimum", J::Number(worst->value)},
                                 {"after", wall}, {"radius", J::Number(radius)},
                                 {"displacement", J::Number(Length(displacement))}});
                            break;
                        }
                    }
                }
                if (groundModelOK) break;
            }
            auto cut  = *worst;
            cut.value = cut.row.Value(state, offset[K()]);
            groundCut(cut);
            ++groundAdded;
            std::fill(solution.begin(), solution.end(), 0.);
            if (restore)
                solution[n] = base.violation + 1e-10;
            active.clear();
        }
        qpOK = qpOK && groundModelOK;
        if (!groundModelOK && quadratic.solved)
            quadratic.failure = "ground_cut_limit";
        if (m_reducedVerify)
        {
            // Replay the small affine problem without repeating field solves.
            // Native little-endian uint64 dimensions followed by f64 arrays:
            // H, g, feasible initial point, then [a,b] for every cut.
            const std::string file = "reduced-model-" + stageTag + "-" + std::to_string(iterations) +
                                     "-" + std::to_string(failedModels) + ".bin";
            std::ofstream out(file, std::ios::binary | std::ios::trunc);
            const std::uint64_t shape[2] = {
                static_cast<std::uint64_t>(variables),
                static_cast<std::uint64_t>(cuts.size())};
            out.write(reinterpret_cast<const char *>(shape), sizeof(shape));
            auto write = [&](const Vector &v) {
                out.write(reinterpret_cast<const char *>(v.data()),
                          v.size() * sizeof(double));
            };
            for (const auto &row : hessian) write(row);
            write(gradient);
            Vector initial(variables);
            if (restore) initial[n] = base.violation + 1e-10;
            write(initial);
            for (const auto &cut : cuts)
            {
                write(cut.a);
                out.write(reinterpret_cast<const char *>(&cut.b), sizeof(double));
            }
            out.close();
            Require(out.good(), "failed to write reduced affine replay packet");
            emit("model_packet", {{"path", J::Quote(file)},
                 {"variables", J::Number(variables)},
                 {"cuts", J::Number(cuts.size())}});
        }
        auto predicted = base.r;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < dimension; ++j)
                predicted[j] += solution[i] * columns[i].r[j];
        double predictedViolation = 0.;
        double unscaledViolation = 0.;
        for (int j = 0; j < base.bounds.size(); ++j)
        {
            const auto &b = base.bounds[j];
            double change = 0.;
            double rawChange = 0.;
            for (int i = 0; i < n; ++i)
            {
                change += solution[i] * columns[i].observations[j];
                if (quadratic.solution.size() == variables)
                    rawChange += quadratic.solution[i] * columns[i].observations[j];
            }
            predictedViolation =
                std::max(predictedViolation,
                         std::max((b.value - b.upper) / b.width + change,
                                  (b.lower - b.value) / b.width - change));
            unscaledViolation = std::max(unscaledViolation,
                std::max((b.value - b.upper)/b.width + rawChange,
                         (b.lower - b.value)/b.width - rawChange));
        }
        auto predictedChange = predicted;
        for (int j = 0; j < dimension; ++j) predictedChange[j] -= base.r[j];
        const double slope = Dot(base.r, mass.Apply(predictedChange)) /
                             (residualScale * residualScale);
        Vector originalMeritGradient(n), actionUncertainty(n), observationUncertainty(n);
        for (int i = 0; i < n; ++i)
        {
            originalMeritGradient[i] = Dot(columns[i].r, dualResidual) /
                                       (residualScale * residualScale);
            actionUncertainty[i] = base.merit * columns[i].error *
                                      std::max(1., dualNorm(columns[i].r));
            observationUncertainty[i] = columns[i].observationError;
        }
        std::vector<std::string> activeRecords;
        std::vector<std::string> basisNames;
        for (const auto &column : columns) basisNames.push_back(J::Quote(column.name));
        for (int id : active)
        {
            const auto &cut = cuts[id];
            std::vector<std::pair<std::string,std::string>> fields{
                {"kind",J::Quote(cut.kind)},{"row",J::Numbers(cut.a)},{"rhs",J::Number(cut.b)}};
            if (cut.observation >= 0)
            {
                const auto &b = base.bounds[cut.observation];
                fields.emplace_back("family",J::Number(b.family));
                fields.emplace_back("point",J::Numbers(b.row.coordinate));
                fields.emplace_back("projected",J::Boolean(b.projected));
                fields.emplace_back("sign",J::Number(cut.sign));
                fields.emplace_back("value",J::Number(b.value));
            }
            else if (groundLocations.count(id))
                fields.emplace_back("point",J::Numbers(groundLocations.at(id).row.coordinate));
            else if (domainLocations.count(id))
            {
                const auto &p = domainLocations.at(id);
                fields.emplace_back("point", J::Numbers(p.coordinate));
                fields.emplace_back("inventory", J::Quote(p.inventory));
                fields.emplace_back("strict", J::Boolean(p.strict));
            }
            activeRecords.push_back(J::Object(fields));
        }
        if (qpOK)
        {
            activeWitnesses.clear();
            std::set<std::pair<int,int>> seen;
            for (int id : active)
                if (cuts[id].observation >= 0)
                {
                    const auto &b = base.bounds[cuts[id].observation];
                    if (seen.emplace(b.family, static_cast<int>(std::floor(4.*(b.row.coordinate[0]-left)))).second)
                        activeWitnesses.push_back(b);
                }
            if (restore && unscaledViolation > .5 * base.violation)
                basisBudget = std::min(32, basisBudget + 3);
        }
        emit("direction",
             {{"restoration", J::Boolean(restore)},
              {"residual_weight", J::Number(residualWeight)},
              {"restoration_balance", J::Number(restorationWeight)},
              {"inward_qualification_search_margin", J::Number(qualificationBuffer)},
              {"qualification_funnel", J::Number(restore ? base.violation * 1.01 + 1e-7 : 0.)},
              {"next_basis_budget", J::Number(basisBudget)},
              {"basis_dimension", J::Number(n)},
              {"basis_names", J::List(basisNames)},
              {"response_rank", J::Number(actionOrthogonal.size())},
              {"qp_solved", J::Boolean(qpOK)},
              {"qp_steps", J::Number(qpSteps)},
              {"redundant_rows", J::Number(quadratic.redundantRows)},
              {"qp_failure", J::Quote(quadratic.failure)},
              {"unscaled_affine_kkt", J::Numbers(Vector{quadratic.primal,
                   quadratic.stationarity, quadratic.complementarity})},
              {"active_constraints", J::Number(active.size())},
              {"active",J::List(activeRecords)},
              {"ground_cuts_added", J::Number(groundAdded)},
              {"ground_model_admissible", J::Boolean(groundModelOK)},
              {"radius", J::Number(radius)},
              {"displacement", J::Number(Length(displacement))},
              {"predicted_norm", J::Number(dualNorm(predicted))},
              {"predicted_violation", J::Number(predictedViolation)},
              {"unscaled_predicted_violation", J::Number(unscaledViolation)},
              {"merit_slope", J::Number(slope)},
              {"action_uncertainty", J::Number(maxActionError)},
              {"original_merit_gradient", J::Numbers(originalMeritGradient)},
              {"two_scale_gradient_error_estimate", J::Numbers(actionUncertainty)},
              {"observation_action_relative_error", J::Numbers(observationUncertainty)}});
        bool accepted = false;
        bool decomposed = false;
        if (qpOK && Length(displacement) > 0. &&
            (restore ? (RANSWallReduced::QualificationDecrease(base.violation, predictedViolation, 1., 1e-6) ||
                        slope < -1e-8 * base.merit * base.merit)
                     : slope < -1e-8 * base.merit * base.merit))
        {
            for (double fraction : {1., .5, .25, .125, .0625, .03125, .015625, .0078125, .00390625})
            {
                auto guess = state;
                Apply(guess, rows, displacement, kScale * fraction);
                if (!admissible(guess))
                    continue;
                auto trial = evaluate(guess, innerTolerance, m_reducedHUpdates,
                                      "reduced_trial");
                if (!trial.valid)
                    continue;
                auto model = base.r;
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < dimension; ++j)
                        model[j] += fraction * solution[i] * columns[i].r[j];
                auto error = trial.r;
                for (int i = 0; i < dimension; ++i)
                    error[i] -= model[i];
                double modelViolation = 0.;
                for (int j = 0; j < base.bounds.size(); ++j)
                {
                    const auto &b = base.bounds[j];
                    double change = 0.;
                    for (int i = 0; i < n; ++i)
                        change +=
                            fraction * solution[i] * columns[i].observations[j];
                    modelViolation = std::max(
                        modelViolation,
                        std::max((b.value - b.upper) / b.width + change,
                                 (b.lower - b.value) / b.width - change));
                }
                bool filterOK = true;
                for (auto [merit, violation] : filter)
                    if (trial.merit >= merit * (1. - 1e-6) &&
                        !RANSWallReduced::QualificationDecrease(violation, trial.violation, 1., 1e-5))
                        filterOK = false;
                const bool progress =
                    restore
                        ? (RANSWallReduced::QualificationDecrease(base.violation, trial.violation, fraction, 1e-3) ||
                           (trial.merit < base.merit * (1. - 1e-3 * fraction) &&
                            trial.violation <= base.violation * 1.01 + 1e-7))
                        : (trial.merit < base.merit && trial.violation < 1e-8);
                const bool take = filterOK && progress;
                if (!decomposed && dualNorm(error) > .2 * std::max(1., base.merit))
                {
                    auto component = [&](const Vector &s) {
                        return extract(ResidualResolutionComponents(s, history, alpha, time, data, false)
                                           .terms.at("inner_k_destruction"));
                    };
                    const auto baseD = component(state);
                    auto predictedD = baseD;
                    for (int i = 0; i < n; ++i)
                    {
                        const auto shiftedD = component(columns[i].shifted.response.state);
                        for (int j = 0; j < dimension; ++j)
                            predictedD[j] += fraction * solution[i] *
                                (shiftedD[j] - baseD[j]) / columns[i].h;
                    }
                    auto destructionError = component(trial.response.state);
                    auto remainder = error;
                    for (int j = 0; j < dimension; ++j)
                    {
                        destructionError[j] -= predictedD[j];
                        remainder[j] -= destructionError[j];
                    }
                    emit("prediction_decomposition",
                         {{"total", J::Number(dualNorm(error))},
                          {"inner_destruction", J::Number(dualNorm(destructionError))},
                          {"remaining", J::Number(dualNorm(remainder))},
                          {"base_b", J::Number(base.domain.minimumB)},
                          {"trial_b", J::Number(trial.domain.minimumB)},
                          {"base_lambda", J::Number(base.domain.minimumLambda)},
                          {"trial_lambda", J::Number(trial.domain.minimumLambda)}});
                    decomposed = true;
                }
                emit("trial",
                     {{"fraction", J::Number(fraction)},
                      {"accepted", J::Boolean(take)},
                      {"norm", J::Number(trial.merit)},
                      {"full_norm", J::Number(Maximum(trial.response.norms))},
                      {"violation", J::Number(trial.violation)},
                      {"families", J::Numbers(trial.family)},
                      {"qualification_inventory", OuterViolationR5(trial.bounds).report},
                      {"predicted_norm", J::Number(dualNorm(model))},
                      {"model_error", J::Number(dualNorm(error))},
                      {"predicted_violation", J::Number(modelViolation)},
                      {"non_k_error", J::Number(trial.response.error)},
                      {"displacement",
                       J::Number(fraction * Length(displacement))},
                      {"physical_k_displacement",
                       J::Number(kScale * fraction * Length(displacement))}});
                if (!take)
                    continue;
                const double oldViolation = base.violation;
                previous = displacement;
                for (double &x : previous)
                    x *= fraction;
                state               = trial.response.state;
                m_lastFeasibleState = state;
                base                = std::move(trial);
                filter.emplace_back(base.merit, base.violation);
                const double predictionError = dualNorm(error);
                const bool agreement = predictionError < .2 * std::max(1., base.merit)
                    && std::abs(base.violation - modelViolation) <
                       .2 * std::max(1., base.violation);
                radius = std::min(1.,
                                  radius * (fraction == 1. && agreement ? 1.8 : .6));
                ++iterations;
                failedModels = 0;
                report(base, "accepted");
                if (restore && base.violation > 0.)
                {
                    const double restored = (oldViolation - base.violation) /
                                            oldViolation;
                    qualificationStalls = restored < .01 ? qualificationStalls + 1 : 0;
                    if (qualificationStalls >= 2)
                    {
                        if (restorationWeight > 0.)
                            restorationWeight = std::max(1e-4, .1 * restorationWeight);
                        basisBudget = std::min(32, basisBudget + 3);
                        qualificationStalls = 0;
                        emit("restoration_balance", {{"residual_multiplier", J::Number(restorationWeight)},
                             {"basis_budget", J::Number(basisBudget)},
                             {"last_relative_qualification_improvement", J::Number(restored)}});
                    }
                }
                else
                {
                    qualificationStalls = 0;
                    restorationWeight = m_reducedRestorationWeight;
                }
                accepted = true;
                break;
            }
        }
        if (!accepted)
        {
            ++failedModels;
            radius *= .5;
            innerTolerance = std::max(minimumInnerAccuracy, innerTolerance * .1);
            emit("model_retry", {{"failed_models", J::Number(failedModels)},
                                 {"radius", J::Number(radius)}});
            Require(failedModels < 4 && radius > 1e-6,
                    "reduced search exhausted local response models without "
                    "acceptable progress");
        }
        innerTolerance = std::max(
            minimumInnerAccuracy, std::min(innerTolerance, base.merit * 1e-3));
    }
    Require(Maximum(base.response.norms) <= tolerance,
            "reduced iteration limit before original full-stage convergence");
    const auto audit =
        CompleteRestrictedAudit(state, history, alpha, time, data, true);
    Require(audit.at("r2_qualified")[0] == 1.,
            "reduced algebraic result failed original qualification");
    return state;
}
} // namespace Nektar
