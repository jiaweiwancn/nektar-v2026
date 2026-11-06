///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallReducedModel.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <map>
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

// These are the existing separate Picard field operators at one equilibrated
// base point. They are preconditioners for small response probes, never a
// replacement for the original H residual. No k row or column is stored here.
struct RANSWallEvolutionPatch::NonKPreconditioner
{
    PicardBlock flow, epsilon;
    int nv;
    WallMatrix divergence, inverseGradient;
    std::unique_ptr<Factor> velocityFactor, pressureFactor, epsilonFactor;

    NonKPreconditioner(PicardBlock f, PicardBlock e, int velocitySize)
        : flow(std::move(f)), epsilon(std::move(e)), nv(velocitySize)
    {
        const int np = flow.rows.size()-nv;
        if (nv <= 0 || np <= 0 || flow.matrix.size() != flow.rows.size())
            throw std::logic_error("invalid response preconditioner flow shape");
        WallMatrix a(nv,Vector(nv)), gradient(nv,Vector(np));
        divergence.assign(np,Vector(nv));
        for (int i=0;i<flow.rows.size();++i)
        {
            if (flow.matrix[i].size()!=flow.rows.size())
                throw std::logic_error("nonsquare response preconditioner flow matrix");
            if (i<nv)
            {
                std::copy_n(flow.matrix[i].begin(),nv,a[i].begin());
                std::copy_n(flow.matrix[i].begin()+nv,np,gradient[i].begin());
            }
            else std::copy_n(flow.matrix[i].begin(),nv,divergence[i-nv].begin());
        }
        velocityFactor=std::make_unique<Factor>(a);
        inverseGradient.assign(np,Vector(nv));
        WallMatrix schur(np,Vector(np));
        for (int j=0;j<np;++j)
        {
            Vector rhs(nv);
            for (int i=0;i<nv;++i) rhs[i]=gradient[i][j];
            inverseGradient[j]=velocityFactor->Apply(std::move(rhs));
            const auto column=Multiply(divergence,inverseGradient[j]);
            for (int i=0;i<np;++i)
                schur[i][j]=column[i]-flow.matrix[nv+i][nv+j];
        }
        pressureFactor=std::make_unique<Factor>(schur);
        epsilonFactor=std::make_unique<Factor>(epsilon.matrix);
    }

    Vector FlowStep(const Vector &residual) const
    {
        const int np=flow.rows.size()-nv;
        Vector rhs(nv);
        for (int i=0;i<nv;++i) rhs[i]=-residual.at(flow.rows[i]);
        auto step=velocityFactor->Apply(std::move(rhs));
        auto pressureRhs=Multiply(divergence,step);
        for (int i=0;i<np;++i) pressureRhs[i]+=residual.at(flow.rows[nv+i]);
        const auto pressure=pressureFactor->Apply(std::move(pressureRhs));
        for (int i=0;i<nv;++i)
            for (int j=0;j<np;++j) step[i]-=inverseGradient[j][i]*pressure[j];
        step.insert(step.end(),pressure.begin(),pressure.end());
        auto defect=Multiply(flow.matrix,step);
        Vector selected(flow.rows.size());
        for (int i=0;i<flow.rows.size();++i)
        {
            selected[i]=residual.at(flow.rows[i]);
            defect[i]+=selected[i];
        }
        Require(Maximum(defect)<=1e-9*std::max(1.,Maximum(selected)),
                "inaccurate cached velocity predictor/correction application");
        return step;
    }

    Vector EpsilonStep(const Vector &residual) const
    {
        Vector rhs;
        for (int row:epsilon.rows) rhs.push_back(-residual.at(row));
        return epsilonFactor->Apply(std::move(rhs));
    }
};

std::unique_ptr<RANSWallEvolutionPatch::NonKPreconditioner>
RANSWallEvolutionPatch::PrepareNonKPreconditioner(
    const Vector &state,const Vector &history,double alpha,double time,
    const RANSWallEvolutionData &data,const std::map<int,double> &external) const
{
    auto flow=AssemblePicardBlock(0,state,history,alpha,time,data,external);
    auto epsilon=AssemblePicardBlock(Epsilon(),state,history,alpha,time,data,external);
    const int nv=std::lower_bound(flow.rows.begin(),flow.rows.end(),DifferentialSize())-
                 flow.rows.begin();
    return std::make_unique<NonKPreconditioner>(std::move(flow),std::move(epsilon),nv);
}

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
    int maximumUpdates, const NonKPreconditioner *preconditioner) const
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
            std::vector<int> flowRows;
            Vector flowStep;
            if (preconditioner)
            {
                flowRows=preconditioner->flow.rows;
                flowStep=preconditioner->FlowStep(result.residual);
                ++result.cachedFlowSolves;
            }
            else
            {
                auto correction=PicardFlowCorrection(current,history,alpha,time,data,external);
                flowRows=std::move(correction.block.rows);
                flowStep=std::move(correction.step);
            }
            ++result.flowSolves;
            const double before =
                Dot(HNorms(result.norms, K()), HNorms(result.norms, K()));
            bool accepted = false;
            std::string trialFailure;
            for (int retry = 0; retry < 9; ++retry)
            {
                const double fraction = std::ldexp(1., -retry);
                auto trial            = current;
                Apply(trial, flowRows, flowStep, fraction);
                try
                {
                    if (!admissible(trial))
                        continue;
                    std::vector<int> epsilonRows;
                    Vector step;
                    if (preconditioner)
                    {
                        // Momentum has just changed. Rebuild the original
                        // residual, including the nonlinear matching trace,
                        // before the separate epsilon correction.
                        const auto residual=Residual(trial,history,alpha,time,data);
                        epsilonRows=preconditioner->epsilon.rows;
                        step=preconditioner->EpsilonStep(residual);
                        ++result.cachedEpsilonSolves;
                    }
                    else
                    {
                        auto epsilon=AssemblePicardBlock(
                            Epsilon(),trial,history,alpha,time,data,external);
                        auto rhs=epsilon.residual;
                        for (auto &v:rhs) v=-v;
                        step=Factor(epsilon.matrix).Apply(std::move(rhs));
                        epsilonRows=std::move(epsilon.rows);
                    }
                    ++result.epsilonSolves;
                    Apply(trial, epsilonRows, step, fraction);
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
    emit("target", {{"physical_history_advanced", "false"},
        {"method", J::Quote(m_picardReduced ? "sequential_reduced_response" : "retained_conditional_control")},
        {"qualification_policy", J::Quote(m_picardReduced
            ? "observe_original_qualification_without_enforcing_it_in_algebraic_diagnostic"
            : "retained_control")}});
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
    int maximumIterations, bool continuationDiagnostic) const
{
    if (continuationDiagnostic && !m_picardContinuation)
        throw std::logic_error("unqualified reduced solves require the explicit no-advance continuation wrapper");
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
    int cachedFlowSolves=0,cachedEpsilonSolves=0;
    int restrictedDestructionChecks=0;
    std::unique_ptr<NonKPreconditioner> nonKCache;
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
                        const std::string &purpose, bool reusePreconditioner=false) {
        Evaluation e;
        const auto *preconditioner=reusePreconditioner ? nonKCache.get() : nullptr;
        while (true)
        {
            e.response=EquilibrateNonK(guess,history,alpha,time,data,accuracy,budget,
                                      preconditioner);
            ++responseSolves;
            flowSolves+=e.response.flowSolves;
            epsilonSolves+=e.response.epsilonSolves;
            cachedFlowSolves+=e.response.cachedFlowSolves;
            cachedEpsilonSolves+=e.response.cachedEpsilonSolves;
            emit("non_k_response",
             {{"purpose", J::Quote(purpose)},
              {"base_field_preconditioner",J::Boolean(preconditioner!=nullptr)},
              {"converged", J::Boolean(e.response.converged)},
              {"required", J::Number(accuracy)},
              {"attained", J::Number(e.response.error)},
              {"updates", J::Number(e.response.updates)},
              {"flow_solves", J::Number(e.response.flowSolves)},
              {"epsilon_solves", J::Number(e.response.epsilonSolves)},
              {"cached_flow_solves",J::Number(e.response.cachedFlowSolves)},
              {"cached_epsilon_solves",J::Number(e.response.cachedEpsilonSolves)},
              {"failure", J::Quote(e.response.failure)}});
            if (!e.response.converged && e.response.state!=guess)
                throw std::logic_error("failed cached non-k solve changed its caller state");
            if (e.response.converged || !preconditioner) break;
            // A stale/unsuitable field preconditioner is a numerical retry,
            // not an H solution. Retry from the exact incoming guess using
            // the retained current-operator separate field solves.
            emit("response_preconditioner_retry",{{"purpose",J::Quote(purpose)},
                {"failure",J::Quote(e.response.failure)}});
            preconditioner=nullptr;
        }
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
        {"maximum_radius", J::Number(m_reducedMaximumRadius)},
        {"response_method", J::Quote("one-sided differences of converged sequential non-k solves")},
        {"half_response_warm_start", J::Quote("interpolate converged non-k coefficients; preserve exact requested k")},
        {"response_field_preconditioner_reuse",J::Boolean(m_reducedResponseCache)},
        {"response_chord_update_budget",J::Number(m_reducedResponseReuseUpdates)},
        {"directional_refresh_budget",J::Number(m_reducedDirectionalRefreshUpdates)},
        {"directional_refresh_policy",J::Quote("optional same-base measured action interpolation before re-solving the residual and individual qualification model; not a derivative certificate")},
        {"response_chord_tight_checks","true"},
        {"response_chord_check_accuracy_limit",J::Number(std::min(tolerance,m_reducedHTolerance))},
        {"response_chord_policy",J::Quote("opt-in historical response secants; fresh two-size signed direction checks; every actual trial converges original H; fresh model after failure, branch change, enlargement or age limit")},
        {"response_preconditioner_policy",J::Quote("probe-only separate base velocity/pressure and epsilon factors; original H convergence required; fresh-operator retry on failure; canonical check and actual trials use current operators")},
        {"initial_basis_budget", J::Number(m_reducedBasis)},
        {"initial_restoration_weight", J::Number(m_reducedRestorationWeight)},
        {"response_krylov_budget", J::Number(m_reducedKrylov)},
        {"active_witness_selection", J::Quote("original_observation_order_per_family_and_region")},
        {"previous_direction_priority", J::Quote("after_conditional_krylov_and_rate_seeds")},
        {"probe_basis", J::Quote("orthogonal field directions with feasible sign/size selection; inward fallback")},
        {"minimum_orthogonal_probe_step", J::Number(1e-10)},
        {"canonical_half_branch_check", J::Quote("always; optional verification cannot change column selection")},
        {"combined_direction_checks", J::Boolean(m_reducedVerify)},
        {"second_order_qualification_correction", J::Boolean(m_reducedSecondOrder)},
        {"constraint_response_secants", J::Boolean(m_reducedConstraintSecant)},
        {"adaptive_qualified_response_space", J::Boolean(m_reducedAdaptiveSpace)},
        {"nonlinear_destruction_ray", J::Boolean(m_reducedSourceRay)},
        {"restricted_destruction_value",J::Boolean(m_reducedFastDestruction)},
        {"nonlinear_destruction_refinements", J::Number(m_reducedSourceRefinements)},
        {"interior_source_model_anchors", J::Boolean(m_reducedInteriorSourceAnchor)},
        {"model_face_solver", J::Quote("Euclidean null space through " +
            std::to_string(RANSWallReduced::DirectFaceThreshold) +
            " model variables; checked augmented orthonormal-face solve above; null-space fallback")},
        {"model_direct_face_threshold", J::Number(RANSWallReduced::DirectFaceThreshold)},
        {"model_cut_filter", J::Quote("large slack-free models: conservative field-Gram ball bound; final scan of every original frozen cut")},
        {"qualification_correction_budget", J::Number(m_reducedConstraintSecant ? 6 : 3)},
        {"continuation_algebraic_diagnostic", J::Boolean(continuationDiagnostic)},
        {"qualification_policy", J::Quote(continuationDiagnostic
            ? "record_only_no_physical_acceptance" : "original_bounds_enforced")},
        {"maximum_basis", J::Number(m_reducedMaximumBasis)}, {"response_accuracy", J::Number(innerTolerance)}});
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
            out << J::Object({{"status", J::Quote(continuationDiagnostic
                    ? "algebraic_continuation_trial_only" : "nonlinear_trial_only")},
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
                                                     epsilonSolves})},
              {"cached_field_applications",J::Numbers(std::array<int,2>{
                  cachedFlowSolves,cachedEpsilonSolves})}});
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
        bool orthogonalProbe = false;
        int probeOrientation = 0, probeSearches = 0;
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
    std::vector<Column> chordColumns;
    int chordAge=0, chordOrigin=-1;
    int directionRefreshes=0;
    int modelPacketSerial=0;
    std::string chordBranches;
    auto chordEligible=[&]() {
        return m_reducedResponseReuseUpdates && base.violation==0. && !failedModels &&
            base.merit>1e-6 && !chordColumns.empty() &&
            chordAge<=m_reducedResponseReuseUpdates && chordColumns.size()>=basisBudget &&
            chordBranches==base.domain.profileBranches;
    };
    auto stageRoot = [&]() {
        if (Maximum(base.response.norms) > tolerance) return false;
        const auto audit = CompleteRestrictedAudit(state, history, alpha, time, data, true);
        return continuationDiagnostic || audit.at("r2_qualified")[0] == 1.;
    };
    while (iterations < maximumIterations)
    {
        // Ordinary stages require original equations and qualification. The
        // explicitly labelled continuation wrapper may inspect an unqualified
        // algebraic root, and always rejects physical/history advancement.
        if (stageRoot())
        {
            report(base, continuationDiagnostic ? "converged_algebraic_diagnostic" : "converged");
            return state;
        }
        if (!continuationDiagnostic && base.violation > 0.)
            innerTolerance = std::max(minimumInnerAccuracy,
                std::min(innerTolerance, 1e-7 * base.violation));
        // A nearly cancelling chosen direction can have a much smaller
        // signal than its individual response columns. Tighten both its
        // reference Y(k) and probes, rather than subtract a coarse-H base
        // from accurate probes. Actual line-search trials keep the common
        // declared inexact tolerance.
        const double chordAccuracy=chordEligible()?std::min(innerTolerance,tolerance):innerTolerance;
        if (base.response.error > chordAccuracy)
        {
            auto refined = evaluate(state, chordAccuracy, m_reducedHUpdates,
                                    "refine_base_accuracy");
            if (!refined.valid && chordAccuracy<innerTolerance)
            {
                chordColumns.clear();chordAge=0;
                emit("response_chord_discarded",{{"reason",J::Quote("tight reference response failed: "+refined.response.failure)},
                    {"scope",J::Quote("accepted state unchanged; rebuild fresh responses at the declared outer accuracy")}});
                continue;
            }
            Require(refined.valid, "reduced base accuracy refinement failed: " +
                                       refined.response.failure);
            if (chordAccuracy<innerTolerance)
            {
                auto difference=refined.r;
                for (int j=0;j<rows.size();++j) difference[j]-=base.r[j];
                emit("response_chord_base_accuracy",{{"required",J::Number(chordAccuracy)},
                    {"trial_accuracy",J::Number(innerTolerance)},
                    {"before_non_k_error",J::Number(base.response.error)},
                    {"attained",J::Number(refined.response.error)},
                    {"original_k_norm_change",J::Number(dualNorm(difference))},
                    {"before_norms",J::Numbers(base.response.norms)},
                    {"after_norms",J::Numbers(refined.response.norms)},
                    {"scope",J::Quote("fixed-k non-k accuracy refinement, no physical/history advance")}});
            }
            base = std::move(refined);
            state = base.response.state;
            m_lastFeasibleState = state;
            if (stageRoot())
            {
                report(base, continuationDiagnostic ? "converged_algebraic_diagnostic" : "converged");
                return state;
            }
        }
        const bool restore  = !continuationDiagnostic && base.violation > 0.;
        const int dimension = rows.size();
        const bool chordReuse=chordEligible();
        nonKCache.reset();
        if (m_reducedResponseCache && !chordReuse)
        {
            try
            {
                nonKCache=PrepareNonKPreconditioner(state,history,alpha,time,data,external);
                emit("response_preconditioner_prepared",{
                    {"base_state_sha256",J::Quote(RANSWallR3::Hash(state))},
                    {"velocity_size",J::Number(nonKCache->nv)},
                    {"pressure_gauge_size",J::Number(nonKCache->flow.rows.size()-nonKCache->nv)},
                    {"epsilon_trace_size",J::Number(nonKCache->epsilon.rows.size())}});
            }
            catch (const std::runtime_error &error)
            {
                emit("response_preconditioner_unavailable",{{"failure",J::Quote(error.what())}});
            }
        }
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
        if (!previous.empty())
        {
            // Retain successful nonlinear field information before the large
            // witness/free-mode pool can exhaust a finite response budget.
            // Its current response still needs the same independent checks.
            add("previous", previous);
            emit("previous_direction_seed", {{"field_norm",J::Number(Length(previous))},
                {"basis_budget",J::Number(basisBudget)}});
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
            // Beyond64, retain the established point/trace seed prefix and
            // spend the larger budget on independent free-field directions.
            if (++witness == std::min(32, std::max(5, (basisBudget+1)/2)))
                break;
        }
        if (m_reducedAdaptiveSpace && basisBudget >= 40)
        {
            // Point and trace seeds can exhaust their span before an expanded
            // response budget is filled. In particular, they need not span
            // interior free modes. Append residual-ranked coefficient modes;
            // the same rank, feasible-probe and two-scale H checks below decide
            // which genuinely add information. This is only a seed pool: no
            // full reduced or all-field matrix is assembled. Modal signs here
            // choose an orientation, never certify polynomial positivity.
            const auto dual = mass.Apply(base.r);
            std::vector<int> modes(dimension);
            std::iota(modes.begin(), modes.end(), 0);
            std::stable_sort(modes.begin(), modes.end(), [&](int a, int b) {
                return std::abs(dual[a]) > std::abs(dual[b]);
            });
            for (int i : modes)
            {
                Vector mode(dimension);
                mode[i] = dual[i] > 0. ? -1. : 1.;
                add("expanded_free_mode_" + std::to_string(rows[i]), std::move(mode));
            }
            emit("free_mode_seed_pool", {{"modes", J::Number(dimension)},
                {"basis_budget", J::Number(basisBudget)},
                {"ranking", J::Quote("absolute free-k mass-inverse residual")}});
        }
        std::vector<Vector> orthogonal;
        std::vector<Column> columns;
        double maxActionError = 0.;
        int generatedKrylov = 0;
        if (chordReuse)
        {
            columns=std::move(chordColumns);
            for (auto &column:columns)
            {
                column.observations.clear();
                for (const auto &b:base.bounds)
                    column.observations.push_back(b.Action(
                        b.projected?column.traceAction:column.stateAction,offset[b.block]));
                maxActionError=std::max(maxActionError,column.error);
            }
            emit("response_chord_reused",{{"dimension",J::Number(columns.size())},
                {"age",J::Number(chordAge)}, {"origin_iteration",J::Number(chordOrigin)},
                {"scope",J::Quote("secant model, not current per-column derivative verification; original direction checks and controlled actual H trials required")}});
        }
        for (int seedId = 0; !chordReuse && seedId < seeds.size(); ++seedId)
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
            // Rank testing alone leaves nearly parallel finite-difference
            // columns and can amplify their errors in a cancelling proposal.
            // At contact the opposite orientation can be inward even when q
            // is outward. Near-tangent polynomials can also require a smaller
            // probe because their minimum moves quadratically. Check both
            // signs before resorting to a common inward component, which can
            // make individually accurate actions almost parallel.
            double probeSize = m_reducedStep;
            for (int refinement = 0; refinement < 20 && probeSize >= 1e-10;
                 ++refinement, probeSize *= .5)
            {
                for (int orientation : {1,-1})
                {
                    auto guess = state;
                    Apply(guess,rows,q,kScale*probeSize*orientation);
                    ++column.probeSearches;
                    if (!admissible(guess)) continue;
                    seed = q;
                    for (double &x : seed) x *= orientation;
                    column.orthogonalProbe = true;
                    column.probeOrientation = orientation;
                    break;
                }
                if (column.orthogonalProbe) break;
            }
            if (!column.orthogonalProbe) probeSize = m_reducedStep;
            column.k    = seed;
            Evaluation shifted;
            double h = probeSize;
            for (int retry = 0; retry < 5; ++retry, h *= .5)
            {
                auto guess = state;
                Apply(guess, rows, seed, kScale * h);
                if (!admissible(guess))
                    continue;
                shifted = evaluate(guess, innerTolerance, m_reducedHUpdates,
                                   "response_" + name, true);
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
                    m_reducedHUpdates, "half_response_" + name, true);
                if (!half.valid)
                    half = evaluate(canonicalGuess, halfAccuracy,
                                    m_reducedHUpdates, "retry_half_response_" + name);
                if (!half.valid)
                    continue;
                // Independently check the accelerated warm start on the first
                // usable column of every model. Failure to recover the
                // same response branch rejects this column, rather than mixing
                // different local Y(k) solutions in a reduced model. Since it
                // can reject a column, this algorithmic check is independent
                // of the optional observation/verification selector.
                if (columns.empty())
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
                      {"orthogonal_probe", J::Boolean(column.orthogonalProbe)},
                      {"probe_orientation", J::Number(column.probeOrientation)},
                      {"probe_searches", J::Number(column.probeSearches)},
                      {"independent_component", J::Number(qLength)},
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
                std::min(m_reducedMaximumBasis, basisBudget + failedModels * 3))
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
        double maximumColumnCorrelation = 0.;
        int orthogonalProbes = 0;
        for (int i = 0; i < n; ++i)
        {
            orthogonalProbes += columns[i].orthogonalProbe;
            for (int j = 0; j < i; ++j)
                maximumColumnCorrelation = std::max(maximumColumnCorrelation,
                                                    std::abs(gram[i][j]));
        }
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
        for (int j = 0; !continuationDiagnostic && j < base.bounds.size(); ++j)
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
        int qpSteps = 0, groundAdded = 0, radiusSolves = 0;
        int augmentedFaceSolves = 0, augmentedFaceFallbacks = 0;
        int maximumBallRedundantCuts = 0;
        int warmFaceSolves=0,warmFaceFallbacks=0;
        double radiusMultiplier = 0., rawDisplacement = 0., radiusComplementarity = 0.;
        bool unconstrainedAvailable = false;
        std::string radiusRefinementFailure;
        WallMatrix solvedHessian = hessian, kColumns;
        for (const auto &column : columns) kColumns.push_back(column.k);
        RANSWallReduced::QuadraticResult quadratic;
        auto solveQP = [&]() {
            auto bounded = RANSWallReduced::SolveRadiusQuadratic<Factor>(
                hessian, gradient, cuts, solution, kColumns, radius, 500);
            if (!bounded.unconstrainedAvailable || bounded.failedPositiveQueries)
                emit("radius_query_retry", {{"purpose",J::Quote("main_direction")},
                    {"initial_failure",J::Quote(bounded.unconstrainedFailure)},
                    {"unconstrained_available",J::Boolean(bounded.unconstrainedAvailable)},
                    {"failed_positive_queries",J::Number(bounded.failedPositiveQueries)},
                    {"solved",J::Boolean(bounded.quadratic.solved)},
                    {"failure",J::Quote(bounded.quadratic.failure)},
                    {"radius",J::Number(radius)}, {"displacement",J::Number(bounded.norm)},
                    {"multiplier",J::Number(bounded.multiplier)},
                    {"primal",J::Number(bounded.quadratic.primal)},
                    {"stationarity",J::Number(bounded.quadratic.stationarity)},
                    {"scope",J::Quote("small linear-cut model only; actual trials still converge H")} });
            qpSteps += bounded.steps;
            radiusSolves += bounded.solves;
            augmentedFaceSolves += bounded.augmentedFaceSolves;
            augmentedFaceFallbacks += bounded.augmentedFaceFallbacks;
            maximumBallRedundantCuts = std::max(maximumBallRedundantCuts,bounded.ballRedundantCuts);
            warmFaceSolves+=bounded.warmFaceSolves;
            warmFaceFallbacks+=bounded.warmFaceFallbacks;
            radiusMultiplier = bounded.multiplier;
            rawDisplacement = bounded.rawNorm;
            unconstrainedAvailable = bounded.unconstrainedAvailable;
            radiusComplementarity = bounded.complementarity;
            radiusRefinementFailure = bounded.refinementFailure;
            solvedHessian = std::move(bounded.hessian);
            quadratic = std::move(bounded.quadratic);
            solution = quadratic.solution;
            active = quadratic.active;
            multipliers = quadratic.multipliers;
            return quadratic.solved;
        };
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
                solution[n] = base.violation + qualificationBuffer + 1e-10;
            active.clear();
        }
        qpOK = qpOK && groundModelOK;
        if (!groundModelOK && quadratic.solved)
            quadratic.failure = "ground_cut_limit";
        if (m_reducedVerify)
        {
            // Replay the small affine problem without repeating field solves.
            // Native little-endian uint64 dimensions followed by f64 arrays:
            // H (including the radius multiplier), g, feasible initial point,
            // then [a,b] for every cut.
            const std::string file = "reduced-model-" + stageTag + "-" + std::to_string(iterations) +
                                     "-" + std::to_string(failedModels) +
                                     (m_reducedDirectionalRefreshUpdates ?
                                      "-refresh"+std::to_string(directionRefreshes)+"-model"+std::to_string(modelPacketSerial++) : "") + ".bin";
            std::ofstream out(file, std::ios::binary | std::ios::trunc);
            const std::uint64_t shape[2] = {
                static_cast<std::uint64_t>(variables),
                static_cast<std::uint64_t>(cuts.size())};
            out.write(reinterpret_cast<const char *>(shape), sizeof(shape));
            auto write = [&](const Vector &v) {
                out.write(reinterpret_cast<const char *>(v.data()),
                          v.size() * sizeof(double));
            };
            for (const auto &row : solvedHessian) write(row);
            write(gradient);
            Vector initial(variables);
            if (restore) initial[n] = base.violation + qualificationBuffer + 1e-10;
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
                 {"radius_multiplier", J::Number(radiusMultiplier)},
                 {"inward_qualification_search_margin", J::Number(qualificationBuffer)},
                 {"cuts", J::Number(cuts.size())}});
        }
        auto predicted = base.r;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < dimension; ++j)
                predicted[j] += solution[i] * columns[i].r[j];
        Vector observationDirection(base.bounds.size());
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
            observationDirection[j] = change;
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
        std::vector<std::string> nextWitnessRecords;
        if (qpOK)
        {
            activeWitnesses.clear();
            std::vector<RANSWallReduced::WitnessKey> candidates;
            for (int id : active)
                if (cuts[id].observation >= 0)
                {
                    const int observation=cuts[id].observation;
                    const auto &b=base.bounds[observation];
                    candidates.push_back({observation,b.family,
                        static_cast<int>(std::floor(4.*(b.row.coordinate[0]-left)))});
                }
            for (int id:RANSWallReduced::SelectWitnesses(std::move(candidates)))
            {
                const auto &b=base.bounds[id];
                activeWitnesses.push_back(b);
                nextWitnessRecords.push_back(J::Object({
                    {"observation",J::Number(id)}, {"family",J::Number(b.family)},
                    {"point",J::Numbers(b.row.coordinate)},
                    {"projected",J::Boolean(b.projected)}}));
            }
            if (restore && unscaledViolation > .5 * base.violation)
                basisBudget = std::min(m_reducedMaximumBasis, basisBudget + 3);
        }
        emit("direction",
             {{"restoration", J::Boolean(restore)},
              {"response_model_kind",J::Quote(chordReuse?"historical_secant":"fresh_two_size_responses")},
              {"response_model_age",J::Number(chordReuse?chordAge:0)},
              {"directional_refreshes",J::Number(directionRefreshes)},
              {"gradient_scope",J::Quote(chordReuse?"secant model, not current first-order certificate":"fresh measured original reduced gradient")},
              {"residual_weight", J::Number(residualWeight)},
              {"restoration_balance", J::Number(restorationWeight)},
              {"inward_qualification_search_margin", J::Number(qualificationBuffer)},
              {"qualification_funnel", J::Number(restore ? base.violation * 1.01 + 1e-7 : 0.)},
              {"next_basis_budget", J::Number(basisBudget)},
              {"basis_dimension", J::Number(n)},
              {"basis_names", J::List(basisNames)},
              {"response_rank", J::Number(actionOrthogonal.size())},
              {"orthogonal_probes", J::Number(orthogonalProbes)},
              {"maximum_column_correlation", J::Number(maximumColumnCorrelation)},
              {"qp_solved", J::Boolean(qpOK)},
              {"qp_steps", J::Number(qpSteps)},
              {"radius_model_solves", J::Number(radiusSolves)},
              {"augmented_face_solves", J::Number(augmentedFaceSolves)},
              {"augmented_face_fallbacks", J::Number(augmentedFaceFallbacks)},
              {"maximum_ball_redundant_cuts", J::Number(maximumBallRedundantCuts)},
              {"warm_face_solves",J::Number(warmFaceSolves)},
              {"warm_face_fallbacks",J::Number(warmFaceFallbacks)},
              {"original_model_cuts", J::Number(cuts.size())},
              {"radius_multiplier", J::Number(radiusMultiplier)},
              {"radius_complementarity", J::Number(radiusComplementarity)},
              {"radius_refinement_failure", J::Quote(radiusRefinementFailure)},
              {"unconstrained_displacement", J::Number(rawDisplacement)},
              {"unconstrained_available", J::Boolean(unconstrainedAvailable)},
              {"redundant_rows", J::Number(quadratic.redundantRows)},
              {"qp_failure", J::Quote(quadratic.failure)},
              {"unscaled_affine_kkt", J::Numbers(Vector{quadratic.primal,
                   quadratic.stationarity, quadratic.complementarity})},
              {"active_constraints", J::Number(active.size())},
              {"active",J::List(activeRecords)},
              {"next_witness_policy",J::Quote("original_observation_order_per_family_and_region")},
              {"next_active_witnesses",J::List(nextWitnessRecords)},
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
        bool directionModelRefreshed=false;
        auto checkChordDirection=[&](const Vector &coefficients,const std::string &purpose) {
            const std::string incomingStateHash=RANSWallR3::Hash(state);
            const std::string incomingHistoryHash=RANSWallR3::Hash(history);
            auto z=coefficients;
            const auto field=RANSWallReduced::BoundCorrection(z,kColumns,
                std::numeric_limits<double>::max());
            if (!(Length(field)>1e-8)) return false;
            const double f=std::min(1e-3,m_reducedStep/Length(field));
            auto guess=state;
            Apply(guess,rows,field,kScale*f);
            auto coarse=evaluate(guess,chordAccuracy,m_reducedHUpdates,"chord_direction");
            guess=state;
            Apply(guess,rows,field,kScale*f*.5);
            auto fine=evaluate(guess,chordAccuracy,m_reducedHUpdates,"half_chord_direction");
            if (!coarse.valid || !fine.valid)
            {
                emit("response_chord_direction_check",{{"purpose",J::Quote(purpose)},
                    {"available","false"},{"passed","false"},
                    {"failure",J::Quote("both controlled H responses are required")}});
                return false;
            }
            Vector action(dimension),fineAction=fine.r,twoSize(dimension),modelError(dimension);
            for (int i=0;i<n;++i)
                for (int j=0;j<dimension;++j) action[j]+=coefficients[i]*columns[i].r[j];
            for (int j=0;j<dimension;++j)
            {
                fineAction[j]=(fineAction[j]-base.r[j])/(.5*f);
                twoSize[j]=(coarse.r[j]-base.r[j])/f-fineAction[j];
                modelError[j]=action[j]-fineAction[j];
            }
            Vector ds=coarse.domain.effectiveState,dt=coarse.domain.traceState,
                   hs=fine.domain.effectiveState,ht=fine.domain.traceState;
            for (int j=0;j<Size();++j)
            {
                ds[j]=(ds[j]-base.domain.effectiveState[j])/f;
                dt[j]=(dt[j]-base.domain.traceState[j])/f;
                hs[j]=(hs[j]-base.domain.effectiveState[j])/(.5*f);
                ht[j]=(ht[j]-base.domain.traceState[j])/(.5*f);
            }
            double observationTwoSize=0.,observationModelError=0.;
            double boundScaledError=0.,boundScaledTwoSize=0.;
            std::map<std::array<int,4>,std::pair<int,double>> witnesses;
            for (int j=0;j<base.bounds.size();++j)
            {
                const auto &b=base.bounds[j];
                const double observed=b.Action(b.projected?ht:hs,offset[b.block]);
                const double coarseObserved=b.Action(b.projected?dt:ds,offset[b.block]);
                double predictedObservation=0.;
                for (int i=0;i<n;++i) predictedObservation+=coefficients[i]*columns[i].observations[j];
                const double scale=std::max(1.,std::abs(observed));
                observationTwoSize=std::max(observationTwoSize,std::abs(coarseObserved-observed)/scale);
                observationModelError=std::max(observationModelError,std::abs(predictedObservation-observed)/scale);
                boundScaledError=std::max(boundScaledError,std::abs(predictedObservation-observed)/b.width);
                boundScaledTwoSize=std::max(boundScaledTwoSize,std::abs(coarseObserved-observed)/b.width);
                for (int sign:{-1,1})
                {
                    const std::array<int,4> key{{b.family,sign,
                        b.row.coordinate[0]<.5*(lower[0]+upper[0])?0:1,b.projected?1:0}};
                    const double excess=sign*(b.value-(sign<0?b.lower:b.upper))/b.width;
                    auto old=witnesses.find(key);
                    if (old==witnesses.end() || excess>old->second.second) witnesses[key]={j,excess};
                }
            }
            std::vector<std::string> signedObservations;
            for (const auto &[key,witness]:witnesses)
            {
                const int id=witness.first,sign=key[1];
                const auto &b=base.bounds[id];
                const double observed=b.Action(b.projected?ht:hs,offset[b.block]);
                const double coarseObserved=b.Action(b.projected?dt:ds,offset[b.block]);
                double predictedObservation=0.;
                for (int i=0;i<n;++i) predictedObservation+=coefficients[i]*columns[i].observations[id];
                signedObservations.push_back(J::Object({{"observation",J::Number(id)},
                    {"family",J::Number(b.family)}, {"sign",J::Number(sign)},
                    {"region",J::Number(key[2])}, {"block",J::Number(b.block)},
                    {"point",J::Numbers(b.row.coordinate)}, {"projected",J::Boolean(b.projected)},
                    {"signed_base",J::Number(witness.second)},
                    {"predicted_change",J::Number(sign*f*predictedObservation)},
                    {"measured_change",J::Number(sign*f*coarseObserved)},
                    {"half_predicted_change",J::Number(sign*.5*f*predictedObservation)},
                    {"half_measured_change",J::Number(sign*.5*f*observed)}}));
            }
            const double scale=std::max(1.,dualNorm(fineAction));
            const double twoSizeError=dualNorm(twoSize)/scale;
            const double modelRelativeError=dualNorm(modelError)/scale;
            const double coarseSlope=(coarse.merit*coarse.merit-base.merit*base.merit)/(2.*f);
            const double halfSlope=(fine.merit*fine.merit-base.merit*base.merit)/f;
            const bool branchChange=coarse.domain.profileBranches!=base.domain.profileBranches ||
                                    fine.domain.profileBranches!=base.domain.profileBranches;
            const bool accurate=!branchChange && twoSizeError<=.03 && observationTwoSize<=.05;
            // In bound-width units, a small raw observation discrepancy can
            // still be large enough to change the active qualification face.
            // A reliable measured action supplies new direction-selection
            // information at this SAME fixed-k base. It does not accept a
            // state or merely repair the constraints after choosing a step.
            const bool modelConsistent=modelRelativeError<=.03 && boundScaledError<=.1;
            const bool passed=accurate &&
                modelRelativeError<=.15 && observationModelError<=.1 && coarseSlope<0. && halfSlope<0. &&
                (!m_reducedDirectionalRefreshUpdates || modelConsistent);
            emit("response_chord_direction_check",{{"purpose",J::Quote(purpose)},
                {"available","true"},{"passed",J::Boolean(passed)},
                {"age",J::Number(chordAge)}, {"required",J::Number(chordAccuracy)},
                {"attained",J::Number(coarse.response.error)},
                {"half_attained",J::Number(fine.response.error)},
                {"two_size_residual_error",J::Number(twoSizeError)},
                {"two_size_observation_error",J::Number(observationTwoSize)},
                {"model_residual_error",J::Number(modelRelativeError)},
                {"model_observation_error",J::Number(observationModelError)},
                {"bound_scaled_observation_model_error",J::Number(boundScaledError)},
                {"bound_scaled_observation_two_size_error",J::Number(boundScaledTwoSize)},
                {"fraction",J::Number(f)}, {"field_norm",J::Number(Length(field))},
                {"measured_action_norm",J::Number(dualNorm(fineAction))},
                {"half_residual_signal",J::Number(.5*f*dualNorm(fineAction))},
                {"coarse_violation",J::Number(coarse.violation)},
                {"half_violation",J::Number(fine.violation)},
                {"signed_observations",J::List(signedObservations)},
                {"predicted_merit_slope",J::Number(Dot(originalMeritGradient,coefficients))},
                {"actual_merit_slope",J::Number(coarseSlope)},
                {"half_merit_slope",J::Number(halfSlope)},
                {"model_consistent_in_bound_units",J::Boolean(modelConsistent)},
                {"directional_refreshes",J::Number(directionRefreshes)},
                {"profile_branch_change",J::Boolean(branchChange)}});
            if (m_reducedDirectionalRefreshUpdates && !modelConsistent && accurate &&
                boundScaledTwoSize<=.05 && directionRefreshes<m_reducedDirectionalRefreshUpdates)
            {
                try
                {
                    RANSWallReduced::Matrix rActions,stateActions,traceActions;
                    for (const auto &column:columns)
                    {
                        rActions.push_back(column.r);
                        stateActions.push_back(column.stateAction);
                        traceActions.push_back(column.traceAction);
                    }
                    const bool updated=RANSWallReduced::UpdateSecantActions(
                            rActions,kColumns,coefficients,fineAction,1e-8) &&
                        RANSWallReduced::UpdateSecantActions(stateActions,kColumns,coefficients,hs,1e-8) &&
                        RANSWallReduced::UpdateSecantActions(traceActions,kColumns,coefficients,ht,1e-8);
                    if (!updated) throw std::runtime_error("insufficient measured direction-action signal");
                    auto defect=[&](const RANSWallReduced::Matrix &actions,Vector observed) {
                        for (int i=0;i<n;++i)
                            for (int j=0;j<observed.size();++j) observed[j]-=coefficients[i]*actions[i][j];
                        return observed;
                    };
                    const double residualDefect=dualNorm(defect(rActions,fineAction))/scale;
                    const double stateDefect=Length(defect(stateActions,hs))/std::max(1.,Length(hs));
                    const double traceDefect=Length(defect(traceActions,ht))/std::max(1.,Length(ht));
                    if (std::max({residualDefect,stateDefect,traceDefect})>1e-8)
                        throw std::runtime_error("inaccurate direction-action interpolation");
                    auto refreshed=columns;
                    for (int i=0;i<n;++i)
                    {
                        auto &column=refreshed[i];
                        column.r=std::move(rActions[i]);
                        column.stateAction=std::move(stateActions[i]);
                        column.traceAction=std::move(traceActions[i]);
                        // The independent k action is exact, regardless of
                        // finite-difference cancellation in the other fields.
                        for (int j=offset[K()];j<offset[Epsilon()];++j)
                            column.stateAction[j]=column.traceAction[j]=0.;
                        for (int j=0;j<dimension;++j)
                            column.stateAction[rows[j]]=column.traceAction[rows[j]]=kScale*column.k[j];
                        column.shifted=Evaluation{};
                        column.observations.clear();
                    }
                    chordColumns=std::move(refreshed);
                    ++directionRefreshes;
                    directionModelRefreshed=true;
                    emit("response_direction_refreshed",{{"purpose",J::Quote(purpose)},
                        {"refresh",J::Number(directionRefreshes)},
                        {"budget",J::Number(m_reducedDirectionalRefreshUpdates)},
                        {"before_state_sha256",J::Quote(incomingStateHash)},
                        {"state_sha256",J::Quote(RANSWallR3::Hash(state))},
                        {"before_history_sha256",J::Quote(incomingHistoryHash)},
                        {"history_sha256",J::Quote(RANSWallR3::Hash(history))},
                        {"residual_interpolation_defect",J::Number(residualDefect)},
                        {"state_interpolation_defect",J::Number(stateDefect)},
                        {"trace_interpolation_defect",J::Number(traceDefect)},
                        {"previous_bound_scaled_error",J::Number(boundScaledError)},
                        {"two_size_bound_scaled_error",J::Number(boundScaledTwoSize)},
                        {"required",J::Number(chordAccuracy)},
                        {"attained",J::Number(coarse.response.error)},
                        {"half_attained",J::Number(fine.response.error)},
                        {"scope",J::Quote("rebuild all signed observation rows and re-solve the constrained model before any actual trial; accepted state, k and physical history unchanged")}});
                }
                catch (const std::runtime_error &error)
                {
                    emit("response_direction_refresh_failed",{{"failure",J::Quote(error.what())}});
                    return false;
                }
                return false;
            }
            return passed;
        };
        if (chordReuse && qpOK && !checkChordDirection(solution,"linear"))
        {
            if (directionModelRefreshed) continue;
            chordColumns.clear();chordAge=0;
            emit("response_chord_discarded",{{"reason",J::Quote("fresh signed direction or accuracy check failed; rebuild responses without changing caller state")}});
            continue;
        }
        if (!chordReuse && m_reducedVerify && qpOK && base.merit > 1e-6 && Length(displacement) > 1e-8)
        {
            // Check the actual combined direction before its larger trials.
            // Individually accurate columns need not predict a cancelling
            // combination accurately, especially at constitutive branch changes.
            const double fraction = std::min(1e-3, m_reducedStep / Length(displacement));
            auto guess = state;
            Apply(guess, rows, displacement, kScale * fraction);
            auto coarse = evaluate(guess, innerTolerance, m_reducedHUpdates,
                                   "combined_direction");
            guess = state;
            Apply(guess, rows, displacement, kScale * fraction * .5);
            auto fine = evaluate(guess, innerTolerance, m_reducedHUpdates,
                                 "half_combined_direction");
            if (coarse.valid && fine.valid)
            {
                auto coarseAction = coarse.r, fineAction = fine.r;
                Vector coarseError(dimension), fineError(dimension);
                for (int j = 0; j < dimension; ++j)
                {
                    coarseAction[j] = (coarseAction[j]-base.r[j])/fraction;
                    fineAction[j] = (fineAction[j]-base.r[j])/(fraction*.5);
                    coarseError[j] = coarseAction[j]-predictedChange[j];
                    fineError[j] = fineAction[j]-predictedChange[j];
                }
                Vector ds = coarse.response.state, dt = coarse.domain.traceState,
                       hs = fine.response.state, ht = fine.domain.traceState;
                for (int j = 0; j < Size(); ++j)
                {
                    ds[j] -= state[j]; hs[j] -= state[j];
                    dt[j] -= base.domain.traceState[j];
                    ht[j] -= base.domain.traceState[j];
                }
                double coarseObservationError = 0., fineObservationError = 0.,
                       coarseChangeError = 0., fineChangeError = 0.;
                std::vector<std::string> observations;
                std::set<std::pair<int,int>> selected;
                for (int id : active)
                    if (cuts[id].observation >= 0)
                        selected.emplace(cuts[id].observation,cuts[id].sign);
                for (int j = 0; j < base.bounds.size(); ++j)
                {
                    const auto &b = base.bounds[j];
                    const double measured = b.Action(b.projected ? dt : ds,offset[b.block]);
                    const double halfMeasured = b.Action(b.projected ? ht : hs,offset[b.block]);
                    const double fullError = measured-fraction*observationDirection[j],
                                 halfError = halfMeasured-fraction*.5*observationDirection[j];
                    coarseChangeError = std::max(coarseChangeError,std::abs(fullError));
                    fineChangeError = std::max(fineChangeError,std::abs(halfError));
                    coarseObservationError = std::max(coarseObservationError,
                        std::abs(fullError/fraction)/std::max(1.,std::abs(measured/fraction)));
                    fineObservationError = std::max(fineObservationError,
                        std::abs(halfError/(fraction*.5))/std::max(1.,std::abs(halfMeasured/(fraction*.5))));
                    for (int sign : {-1,1}) if (selected.count({j,sign}))
                        observations.push_back(J::Object({{"family",J::Number(b.family)},
                            {"sign",J::Number(sign)},{"point",J::Numbers(b.row.coordinate)},
                            {"projected",J::Boolean(b.projected)},
                            {"predicted_change",J::Number(sign*fraction*observationDirection[j])},
                            {"measured_change",J::Number(sign*measured)},
                            {"half_predicted_change",J::Number(sign*fraction*.5*observationDirection[j])},
                            {"half_measured_change",J::Number(sign*halfMeasured)}}));
                }
                emit("combined_direction_check", {{"available","true"},
                    {"role",J::Quote(restore ? "restoration" : "residual_reduction")},
                    {"fraction",J::Number(fraction)},
                    {"required",J::Number(innerTolerance)},
                    {"attained",J::Number(coarse.response.error)},
                    {"half_attained",J::Number(fine.response.error)},
                    {"relative_residual_error",J::Number(dualNorm(coarseError)/std::max(1.,dualNorm(coarseAction)))},
                    {"half_relative_residual_error",J::Number(dualNorm(fineError)/std::max(1.,dualNorm(fineAction)))},
                    {"observation_relative_error",J::Number(coarseObservationError)},
                    {"half_observation_relative_error",J::Number(fineObservationError)},
                    {"observation_change_error",J::Number(coarseChangeError)},
                    {"half_observation_change_error",J::Number(fineChangeError)},
                    {"predicted_merit_slope",J::Number(slope)},
                    {"actual_merit_slope",J::Number((coarse.merit*coarse.merit-base.merit*base.merit)/(2.*fraction))},
                    {"half_merit_slope",J::Number((fine.merit*fine.merit-base.merit*base.merit)/fraction)},
                    {"full_inventory_violation",J::Number(coarse.violation)},
                    {"half_full_inventory_violation",J::Number(fine.violation)},
                    {"profile_branch_change",J::Boolean(coarse.domain.profileBranches!=base.domain.profileBranches)},
                    {"half_profile_branch_change",J::Boolean(fine.domain.profileBranches!=base.domain.profileBranches)},
                    {"signed_active_observations",J::List(observations)}});
            }
            else emit("combined_direction_check", {{"available","false"},
                {"coarse_valid",J::Boolean(coarse.valid)}, {"half_valid",J::Boolean(fine.valid)},
                {"failure",J::Quote("combined probes did not both attain valid controlled H responses")}});
        }
        bool accepted = false;
        bool decomposed = false;
        int correctionEvaluations = 0;
        auto component = [&](const Vector &s) {
            if (!m_reducedFastDestruction)
                return extract(ResidualResolutionComponents(s, history, alpha, time, data, false)
                                   .terms.at("inner_k_destruction"));
            const auto fastStarted=std::chrono::steady_clock::now();
            auto value=RestrictedInnerDestructionValue(s,time,data,rows);
            const double fastSeconds=std::chrono::duration<double>(
                std::chrono::steady_clock::now()-fastStarted).count();
            // This functional guard is independent of the optional verifier.
            // Cover the base, finite affine trial and small tangent sample.
            if (restrictedDestructionChecks<3)
            {
                const auto fullStarted=std::chrono::steady_clock::now();
                const auto full=extract(ResidualResolutionComponents(s,history,alpha,time,data,false)
                    .terms.at("inner_k_destruction"));
                const double fullSeconds=std::chrono::duration<double>(
                    std::chrono::steady_clock::now()-fullStarted).count();
                Require(RANSWallR3::Hash(value)==RANSWallR3::Hash(full),
                    "restricted destruction value differs from the full diagnostic component");
                emit("restricted_destruction_value_check",{{"sample",J::Number(restrictedDestructionChecks++)},
                    {"rows",J::Number(value.size())}, {"bitwise_equal","true"},
                    {"mass_dual_norm",J::Number(dualNorm(value))},
                    {"restricted_seconds",J::Number(fastSeconds)},
                    {"diagnostic_seconds",J::Number(fullSeconds)},
                    {"scope",J::Quote("shared spatial destruction value only; original equations and physical/history acceptance unchanged")}});
            }
            return value;
        };
        Vector sourceBase;
        int sourceModelEvaluations = 0;
        // This is an internal surrogate on the measured affine y response,
        // never an unequilibrated state claimed to be Y(k). Every actual
        // fraction below must still converge H and scan the original fields.
        auto sourceModel = [&](const Vector &coefficients, Vector &model) {
            if (coefficients.size()!=n) throw std::logic_error("source-ray coefficient shape");
            ++sourceModelEvaluations;
            try
            {
                Vector ds(Size()), dk(dimension);
                model=base.r;
                for (int i=0;i<n;++i)
                {
                    for (int j=0;j<Size();++j) ds[j]+=coefficients[i]*columns[i].stateAction[j];
                    for (int j=0;j<dimension;++j)
                    {
                        dk[j]+=coefficients[i]*columns[i].k[j];
                        model[j]+=coefficients[i]*columns[i].r[j];
                    }
                }
                const double length=Length(dk);
                if (!(length>0.)) return true;
                for (int j=offset[K()];j<offset[Epsilon()];++j) ds[j]=0.;
                for (int j=0;j<dimension;++j) ds[rows[j]]=kScale*dk[j];
                const double h=std::min(1e-3,m_reducedStep/length);
                auto predictedState=state, probe=state;
                for (int j=0;j<Size();++j)
                { predictedState[j]+=ds[j]; probe[j]+=h*ds[j]; }
                if (!admissible(predictedState) || !admissible(probe)) return false;
                if (sourceBase.empty()) sourceBase=component(state);
                const auto value=component(predictedState), tangent=component(probe);
                for (int j=0;j<dimension;++j)
                    model[j]+=value[j]-sourceBase[j]-(tangent[j]-sourceBase[j])/h;
                return true;
            }
            catch (const std::runtime_error &) { return false; }
            catch (const std::invalid_argument &) { return false; }
        };
        if (qpOK && Length(displacement) > 0. &&
            (restore ? (RANSWallReduced::QualificationDecrease(base.violation, predictedViolation, 1., 1e-6) ||
                        slope < -1e-8 * base.merit * base.merit)
                     : slope < -1e-8 * base.merit * base.merit))
        {
            struct Proposal { Vector coefficients; double fraction; std::string origin; };
            std::vector<Proposal> proposals;
            auto addProposal=[&](const Vector &z,double f,const std::string &origin) {
                if (!(f>0. && f<=1.)) return;
                auto scaled=z;
                for (double &x:scaled) x*=f;
                for (const auto &old:proposals)
                {
                    auto difference=scaled;
                    for (int i=0;i<n;++i) difference[i]-=old.coefficients[i];
                    if (Length(difference)<1e-14) return;
                }
                proposals.push_back({std::move(scaled),f,origin});
            };
            auto addFraction=[&](double f) {
                addProposal(solution,f,"linear_reduced_direction");
            };
            if (m_reducedSourceRay && !restore && base.merit>1e-6)
            {
                const auto choice=RANSWallReduced::SelectSourceRay(base.merit,
                    [&](double f,double &value) {
                        Vector z(n),model;
                        for (int i=0;i<n;++i) z[i]=f*solution[i];
                        if (!sourceModel(z,model)) return false;
                        value=dualNorm(model); return true;
                    });
                std::vector<std::string> samples;
                for (auto [f,value]:choice.samples)
                    samples.push_back(J::Numbers(Vector{f,value}));
                emit("source_ray_selection", {{"fraction",J::Number(choice.fraction)},
                    {"predicted_norm",J::Number(choice.merit)}, {"base_norm",J::Number(base.merit)},
                    {"evaluations",J::Number(choice.evaluations)}, {"invalid",J::Number(choice.invalid)},
                    {"samples",J::List(samples)},
                    {"scope",J::Quote("surrogate only: nonlinear destruction on affine converged-field responses; actual trials solve H")}});
                if (choice.fraction>0. && m_reducedSourceRefinements)
                {
                    Vector z(n),model;
                    for (int i=0;i<n;++i) z[i]=choice.fraction*solution[i];
                    auto objective=[&](const Vector &q,const Vector &r) {
                        auto coefficients=q;
                        const auto field=RANSWallReduced::BoundCorrection(coefficients,kColumns,
                            std::numeric_limits<double>::max());
                        return dualNorm(r)*dualNorm(r)+1e-3*Dot(field,field)+1e-10*Dot(q,q);
                    };
                    bool changed=false;
                    for (int refinement=0;refinement<m_reducedSourceRefinements;++refinement)
                    {
                        if (!sourceModel(z,model)) break;
                        WallMatrix nextH(n,Vector(n));
                        const int evaluationsBefore=sourceModelEvaluations;
                        const int evaluationBudget=4*n+64;
                        auto boundedModel=[&](const Vector &q,Vector &value) {
                            if (sourceModelEvaluations-evaluationsBefore>=evaluationBudget)
                                return false;
                            return sourceModel(q,value);
                        };
                        auto relativeError=[&](const RANSWallReduced::ModelDerivative &d) {
                            auto difference=d.coarse;
                            for (int j=0;j<dimension;++j) difference[j]-=d.fine[j];
                            return dualNorm(difference)/std::max(1.,dualNorm(d.fine));
                        };
                        auto linearization=RANSWallReduced::LinearizeModelOnRay(
                            z,model,2e-5,m_reducedInteriorSourceAnchor?2:0,.03,
                            boundedModel,relativeError);
                        const double error=linearization.error;
                        const int available=linearization.availableDerivatives;
                        if (linearization.derivativeRefinements)
                        {
                            std::vector<std::string> stencils;
                            for (const auto &s:linearization.stencils)
                                stencils.push_back(J::Object({{"column",J::Number(s.column)},
                                    {"name",J::Quote(columns.at(s.column).name)},
                                    {"halvings",J::Number(s.halvings)},
                                    {"initial_error",J::Number(s.initialError)},
                                    {"final_error",J::Number(s.finalError)},
                                    {"coarse_step",J::Number(s.coarseStep)},
                                    {"fine_step",J::Number(.5*s.coarseStep)}}));
                            emit("source_derivative_precision",{
                                {"refinement",J::Number(refinement)},
                                {"available",J::Boolean(linearization.available)},
                                {"dimension",J::Number(n)},
                                {"refined_columns",J::Number(linearization.refinedDerivatives)},
                                {"stencil_halvings",J::Number(linearization.derivativeRefinements)},
                                {"minimum_coarse_step",J::Number(linearization.minimumStep)},
                                {"minimum_fine_step",J::Number(.5*linearization.minimumStep)},
                                {"initial_error",J::Number(linearization.initialError)},
                                {"final_error",J::Number(error)},
                                {"model_evaluations",J::Number(sourceModelEvaluations-evaluationsBefore)},
                                {"model_evaluation_budget",J::Number(evaluationBudget)},
                                {"columns",J::List(stencils)},
                                {"scope",J::Quote("surrogate derivatives only; accuracy threshold unchanged; actual trials converge H")}});
                        }
                        if (m_reducedInteriorSourceAnchor)
                            emit("source_linearization_anchor",{
                                {"refinement",J::Number(refinement)},
                                {"available",J::Boolean(linearization.available)},
                                {"anchor_reductions",J::Number(linearization.anchorReductions)},
                                {"attempts",J::Number(linearization.attempts)},
                                {"anchor_fraction",J::Number(std::ldexp(1.,-linearization.anchorReductions))},
                                {"available_derivatives",J::Number(available)},
                                {"dimension",J::Number(n)}, {"two_scale_error",J::Number(error)},
                                {"model_evaluations",J::Number(sourceModelEvaluations-evaluationsBefore)},
                                {"model_evaluation_budget",J::Number(evaluationBudget)},
                                {"scope",J::Quote("surrogate-only anchor; incumbent retained; actual trials converge H")}});
                        if (!linearization.available)
                        {
                            emit("source_direction_refinement", {{"refinement",J::Number(refinement)},
                                {"improved","false"},{"failure",J::Quote("unresolved surrogate derivative")},
                                {"available_derivatives",J::Number(available)},
                                {"dimension",J::Number(n)},{"two_scale_error",J::Number(error)}});
                            break;
                        }
                        auto derivative=std::move(linearization.derivatives);
                        auto constant=linearization.value;
                        for (int i=0;i<n;++i)
                            for (int j=0;j<dimension;++j)
                                constant[j]-=linearization.at[i]*derivative[i][j];
                        const auto dualConstant=mass.Apply(constant);
                        Vector nextG(n);
                        for (int i=0;i<n;++i)
                        {
                            const auto dual=mass.Apply(derivative[i]);
                            nextG[i]=Dot(derivative[i],dualConstant)/(residualScale*residualScale);
                            for (int j=0;j<n;++j)
                                nextH[i][j]=Dot(derivative[j],dual)/(residualScale*residualScale)
                                    +1e-3*gram[i][j]+(i==j?1e-10:0.);
                        }
                        RANSWallReduced::RadiusResult refined;
                        try
                        {
                            refined=RANSWallReduced::SolveRadiusQuadratic<Factor>(
                                nextH,nextG,cuts,z,kColumns,radius,500);
                            if (!refined.unconstrainedAvailable || refined.failedPositiveQueries)
                                emit("radius_query_retry", {{"purpose",J::Quote("source_refinement")},
                                    {"initial_failure",J::Quote(refined.unconstrainedFailure)},
                                    {"unconstrained_available",J::Boolean(refined.unconstrainedAvailable)},
                                    {"failed_positive_queries",J::Number(refined.failedPositiveQueries)},
                                    {"solved",J::Boolean(refined.quadratic.solved)},
                                    {"failure",J::Quote(refined.quadratic.failure)},
                                    {"radius",J::Number(radius)}, {"displacement",J::Number(refined.norm)},
                                    {"multiplier",J::Number(refined.multiplier)},
                                    {"primal",J::Number(refined.quadratic.primal)},
                                    {"stationarity",J::Number(refined.quadratic.stationarity)},
                                    {"scope",J::Quote("source surrogate model only; actual trials still converge H")} });
                        }
                        catch (const std::runtime_error &e)
                        {
                            emit("source_direction_refinement", {{"refinement",J::Number(refinement)},
                                {"improved","false"},{"failure",J::Quote(e.what())}});
                            break;
                        }
                        bool improved=false;
                        double chosenFraction=0.,chosenNorm=dualNorm(model),originalSlope=0.;
                        Vector next=z;
                        if (refined.quadratic.solved)
                            for (double f:{1.,.5,.25,.125,.0625,.03125,.015625})
                            {
                                auto q=z;
                                for (int i=0;i<n;++i) q[i]+=f*(refined.quadratic.solution[i]-z[i]);
                                const auto field=RANSWallReduced::BoundCorrection(q,kColumns,
                                    std::numeric_limits<double>::max());
                                if (Length(field)>radius*(1.+1e-12) ||
                                    std::any_of(cuts.begin(),cuts.end(),[&](const Cut &c) {
                                        return Dot(c.a,q)>c.b+1e-8; })) continue;
                                const double trueSlope=Dot(originalMeritGradient,q);
                                if (trueSlope>=-1e-8*base.merit*base.merit) continue;
                                Vector value;
                                if (!sourceModel(q,value)) continue;
                                if (m_reducedInteriorSourceAnchor &&
                                    dualNorm(value)>=dualNorm(model)) continue;
                                if (objective(q,value)>=objective(z,model)*(1.-1e-6*f)) continue;
                                next=std::move(q); chosenFraction=f; chosenNorm=dualNorm(value);
                                originalSlope=trueSlope; improved=true; break;
                            }
                        emit("source_direction_refinement", {{"refinement",J::Number(refinement)},
                            {"improved",J::Boolean(improved)},
                            {"failure",J::Quote(refined.quadratic.failure)},
                            {"dimension",J::Number(n)},{"available_derivatives",J::Number(available)},
                            {"two_scale_error",J::Number(error)},
                            {"augmented_face_solves",J::Number(refined.augmentedFaceSolves)},
                            {"augmented_face_fallbacks",J::Number(refined.augmentedFaceFallbacks)},
                            {"warm_face_solves",J::Number(refined.warmFaceSolves)},
                            {"warm_face_fallbacks",J::Number(refined.warmFaceFallbacks)},
                            {"ball_redundant_cuts",J::Number(refined.ballRedundantCuts)},
                            {"before_norm",J::Number(dualNorm(model))},
                            {"predicted_norm",J::Number(chosenNorm)},
                            {"fraction",J::Number(chosenFraction)},
                            {"original_merit_slope",J::Number(originalSlope)},
                            {"radius",J::Number(radius)},
                            {"scope",J::Quote("constrained nonlinear destruction surrogate in measured response space; not Y(k)")}});
                        if (!improved) break;
                        z=std::move(next); changed=true;
                    }
                    if (changed)
                    {
                        auto q=z;
                        const auto field=RANSWallReduced::BoundCorrection(q,kColumns,
                            std::numeric_limits<double>::max());
                        emit("source_refined_direction", {{"displacement",J::Number(Length(field))},
                            {"radius",J::Number(radius)},
                            {"original_merit_slope",J::Number(Dot(originalMeritGradient,z))},
                            {"coefficients",J::Numbers(z)},
                            {"source_model_evaluations",J::Number(sourceModelEvaluations)}});
                        if (m_reducedVerify && !chordReuse)
                        {
                            const double f=std::min(1e-3,m_reducedStep/Length(field));
                            auto guess=state;
                            Apply(guess,rows,field,kScale*f);
                            auto coarse=evaluate(guess,innerTolerance,m_reducedHUpdates,"source_refined_direction");
                            guess=state;
                            Apply(guess,rows,field,kScale*f*.5);
                            auto fine=evaluate(guess,innerTolerance,m_reducedHUpdates,"half_source_refined_direction");
                            if (coarse.valid && fine.valid)
                            {
                                auto difference=fine.r;
                                Vector action(dimension);
                                for (int i=0;i<n;++i)
                                    for (int j=0;j<dimension;++j) action[j]+=z[i]*columns[i].r[j];
                                for (int j=0;j<dimension;++j)
                                    difference[j]=(difference[j]-base.r[j])/(.5*f)-action[j];
                                emit("source_refined_direction_check", {{"available","true"},
                                    {"required",J::Number(innerTolerance)},
                                    {"attained",J::Number(coarse.response.error)},
                                    {"half_attained",J::Number(fine.response.error)},
                                    {"half_relative_residual_error",J::Number(dualNorm(difference)/std::max(1.,dualNorm(action)))},
                                    {"predicted_merit_slope",J::Number(Dot(originalMeritGradient,z))},
                                    {"actual_merit_slope",J::Number((coarse.merit*coarse.merit-base.merit*base.merit)/(2.*f))},
                                    {"half_merit_slope",J::Number((fine.merit*fine.merit-base.merit*base.merit)/f)}});
                            }
                            else emit("source_refined_direction_check",{{"available","false"}});
                        }
                        if (!chordReuse || checkChordDirection(z,"source_refined"))
                            for (double f:{1.,.5,.25}) addProposal(z,f,"nonlinear_source_refined_direction");
                    }
                }
                if (choice.fraction>0.)
                    for (double f:{choice.fraction,.5*choice.fraction,.25*choice.fraction}) addFraction(f);
            }
            if (directionModelRefreshed) continue;
            for (double f:{1.,.5,.25,.125,.0625,.03125,.015625,.0078125,.00390625}) addFraction(f);
            int proposalIndex=0;
            for (const auto &proposal : proposals)
            {
                ++proposalIndex;
                const double fraction=proposal.fraction;
                auto trialCoefficients=proposal.coefficients;
                auto trialDisplacement=RANSWallReduced::BoundCorrection(trialCoefficients,kColumns,
                    std::numeric_limits<double>::max());
                if (proposal.origin=="linear_reduced_direction")
                {
                    // Preserve the established control's floating-point path.
                    trialDisplacement=displacement;
                    for (double &x:trialDisplacement) x*=fraction;
                }
                auto guess = state;
                Apply(guess, rows, trialDisplacement, kScale);
                if (!admissible(guess))
                    continue;
                auto trial = evaluate(guess, innerTolerance, m_reducedHUpdates,
                                      "reduced_trial");
                if (!trial.valid)
                    continue;
                const double originalTrialLength = Length(trialDisplacement);
                int usedCorrections = 0;
                WallMatrix correctionState, correctionTrace;
                for (const auto &column : columns)
                {
                    correctionState.push_back(column.stateAction);
                    correctionTrace.push_back(column.traceAction);
                }
                auto exactKActions = [&]() {
                    for (int i = 0; i < n; ++i)
                    {
                        for (int j = offset[K()]; j < offset[Epsilon()]; ++j)
                            correctionState[i][j] = correctionTrace[i][j] = 0.;
                        for (int j = 0; j < dimension; ++j)
                            correctionState[i][rows[j]] = correctionTrace[i][rows[j]] =
                                kScale*columns[i].k[j];
                    }
                };
                if (m_reducedConstraintSecant) exactKActions();
                // A tangent reduced step can leave curved qualification
                // surfaces even when its checked first derivative is accurate.
                // Correct the actual trial observations in the same measured
                // response space. This changes k only; each corrected trial
                // must again solve H and scan the complete actual inventory.
                while (m_reducedSecondOrder && !continuationDiagnostic && !restore &&
                       fraction >= (m_reducedSourceRay ? 1e-5 : .25) &&
                       trial.violation > 0. && trial.violation < .1 &&
                       trial.merit < base.merit &&
                       correctionEvaluations < (m_reducedConstraintSecant ? 6 : 3) &&
                       originalTrialLength > 1e-8)
                {
                    const int count = n + 1; // Last variable is an elastic slack.
                    WallMatrix correctionH(count, Vector(count));
                    Vector correctionG(count), initial(count);
                    for (int i = 0; i < n; ++i)
                        for (int j = 0; j < n; ++j)
                            correctionH[i][j] = gram[i][j] + (i == j ? 1e-10 : 0.);
                    correctionH[n][n] = 1e8;
                    const double margin = std::min(1e-6, 1e-4 * base.merit);
                    initial[n] = trial.violation + margin + 1e-10;
                    std::vector<Cut> correctionCuts;
                    auto cut = [&](Vector a, double b, const std::string &kind) {
                        const double norm = Length(a);
                        if (norm < 1e-15) return;
                        for (double &x : a) x /= norm;
                        correctionCuts.push_back({std::move(a), b/norm, -1, 0, kind});
                    };
                    for (const auto &b : trial.bounds)
                        for (int sign : {-1, 1})
                        {
                            Vector a(count);
                            for (int i = 0; i < n; ++i)
                                a[i] = sign * b.Action(b.projected ? correctionTrace[i]
                                                                  : correctionState[i],
                                                       offset[b.block]);
                            a[n] = -1.;
                            cut(std::move(a), -sign * (b.value -
                                (sign < 0 ? b.lower : b.upper))/b.width - margin,
                                "qualification");
                        }
                    for (const auto &p : WallTraceR4(trial.response.state))
                    {
                        Vector a(count);
                        for (int i = 0; i < n; ++i)
                            for (int j = 0; j < p.row.ids.size(); ++j)
                            {
                                const int c = index[p.row.ids[j]];
                                if (c >= 0) a[i] -= kScale * p.row.phi[j] * columns[i].k[c];
                            }
                        cut(std::move(a), p.value - p.error, "continuous_ground");
                    }
                    for (const auto &p : trial.domain.observations)
                    {
                        const auto row = Point(p.block, p.coordinate, p.side);
                        Vector a(count);
                        for (int i = 0; i < n; ++i)
                            a[i] = -row.Value(p.dependentTrace ? correctionTrace[i]
                                                              : correctionState[i],
                                             offset[p.block])/p.scale;
                        cut(std::move(a), p.value/p.scale, "evaluation_domain");
                    }
                    Vector nonnegative(count); nonnegative[n] = -1.;
                    cut(std::move(nonnegative), 0., "elastic_slack");
                    const auto correction = RANSWallReduced::SolveRadiusQuadratic<Factor>(
                        correctionH, correctionG, correctionCuts, initial, kColumns,
                        .5*originalTrialLength, 500);
                    if (!correction.quadratic.solved)
                    {
                        emit("second_order_model_failed", {{"fraction", J::Number(fraction)},
                            {"failure", J::Quote(correction.quadratic.failure)}});
                        break;
                    }
                    auto delta = correction.quadratic.solution;
                    const auto change = RANSWallReduced::BoundCorrection(delta, kColumns,
                        std::numeric_limits<double>::max());
                    auto total = trialDisplacement;
                    for (int i = 0; i < dimension; ++i) total[i] += change[i];
                    if (Length(total) > radius ||
                        Length(change) > .5*originalTrialLength*(1.+1e-12))
                    {
                        emit("second_order_model_failed", {{"fraction", J::Number(fraction)},
                            {"failure", J::Quote("corrected actual field exceeds its displacement bound")}});
                        break;
                    }
                    // Warm start on the same converged H branch, while
                    // prescribing exactly the corrected independent k.
                    auto correctedGuess = trial.response.state;
                    Apply(correctedGuess, rows, change, kScale);
                    if (!admissible(correctedGuess))
                    {
                        emit("second_order_model_failed", {{"fraction", J::Number(fraction)},
                            {"failure", J::Quote("corrected k fails actual ground or evaluation-domain scan")}});
                        break;
                    }
                    ++correctionEvaluations;
                    auto corrected = evaluate(correctedGuess, innerTolerance, m_reducedHUpdates,
                                               "second_order_trial");
                    double predictedQualification = 0.;
                    for (const auto &b : trial.bounds)
                    {
                        double change = 0.;
                        for (int i = 0; i < n; ++i)
                            change += delta[i] * b.Action(b.projected ? correctionTrace[i]
                                                                    : correctionState[i],
                                                          offset[b.block]);
                        predictedQualification = std::max(predictedQualification,
                            std::max((b.value-b.upper)/b.width+change,
                                     (b.lower-b.value)/b.width-change));
                    }
                    const bool improved = corrected.valid && corrected.merit < base.merit &&
                                          corrected.violation < trial.violation;
                    emit("second_order_trial", {{"fraction", J::Number(fraction)},
                        {"valid", J::Boolean(corrected.valid)}, {"improved", J::Boolean(improved)},
                        {"before_norm", J::Number(trial.merit)},
                        {"norm", J::Number(corrected.merit)},
                        {"before_violation", J::Number(trial.violation)},
                        {"violation", J::Number(corrected.violation)},
                        {"predicted_violation", J::Number(predictedQualification)},
                        {"elastic_slack", J::Number(delta[n])},
                        {"inventory_size", J::Number(trial.bounds.size())},
                        {"required", J::Number(innerTolerance)},
                        {"non_k_error", J::Number(corrected.response.error)},
                        {"correction", J::Number(Length(change))},
                        {"correction_radius", J::Number(.5*originalTrialLength)},
                        {"displacement", J::Number(Length(total))}, {"radius", J::Number(radius)}});
                    if (!improved) break;
                    if (m_reducedConstraintSecant)
                    {
                        Vector observedState = corrected.domain.effectiveState,
                               observedTrace = corrected.domain.traceState;
                        for (int j = 0; j < Size(); ++j)
                        {
                            observedState[j] -= trial.domain.effectiveState[j];
                            observedTrace[j] -= trial.domain.traceState[j];
                        }
                        double observationError = 0.;
                        for (const auto &b : trial.bounds)
                        {
                            double predicted = 0.;
                            for (int i = 0; i < n; ++i)
                                predicted += delta[i] * b.Action(b.projected ? correctionTrace[i]
                                                                            : correctionState[i],
                                                                 offset[b.block]);
                            const double measured = b.Action(b.projected ? observedTrace : observedState,
                                                              offset[b.block]);
                            observationError = std::max(observationError,std::abs(measured-predicted));
                        }
                        bool updated = false;
                        std::string failure;
                        try
                        {
                            auto nextState = correctionState, nextTrace = correctionTrace;
                            const Vector coefficients(delta.begin(), delta.begin()+n);
                            updated = RANSWallReduced::UpdateSecantActions(
                                nextState,kColumns,coefficients,observedState,1e-6);
                            const bool traceUpdated = RANSWallReduced::UpdateSecantActions(
                                nextTrace,kColumns,coefficients,observedTrace,1e-6);
                            if (updated != traceUpdated)
                                throw std::logic_error("constraint secant signal tests disagree");
                            correctionState = std::move(nextState);
                            correctionTrace = std::move(nextTrace);
                            // Independent k has a known exact response. Refine
                            // only the other fields and dependent trace target.
                            exactKActions();
                        }
                        catch (const std::runtime_error &error)
                        {
                            updated = false;
                            failure = error.what();
                        }
                        emit("constraint_response_secant", {{"updated",J::Boolean(updated)},
                            {"failure",J::Quote(failure)},
                            {"minimum_field_signal",J::Number(1e-6)},
                            {"field_signal",J::Number(Length(change))},
                            {"maximum_signed_observation_prediction_error",J::Number(observationError)},
                            {"before_violation",J::Number(trial.violation)},
                            {"after_violation",J::Number(corrected.violation)}});
                    }
                    trial = std::move(corrected);
                    trialDisplacement = std::move(total);
                    for (int i = 0; i < n; ++i) trialCoefficients[i] += delta[i];
                    ++usedCorrections;
                }
                auto model = base.r;
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < dimension; ++j)
                        model[j] += trialCoefficients[i] * columns[i].r[j];
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
                            trialCoefficients[i] * columns[i].observations[j];
                    modelViolation = std::max(
                        modelViolation,
                        std::max((b.value - b.upper) / b.width + change,
                                 (b.lower - b.value) / b.width - change));
                }
                bool filterOK = true;
                for (auto [merit, violation] : filter)
                    if (trial.merit >= merit * (1. - 1e-6) &&
                        (continuationDiagnostic ||
                         !RANSWallReduced::QualificationDecrease(violation, trial.violation, 1., 1e-5)))
                        filterOK = false;
                const bool progress =
                    restore
                        ? (RANSWallReduced::QualificationDecrease(base.violation, trial.violation, fraction, 1e-3) ||
                           (trial.merit < base.merit * (1. - 1e-3 * fraction) &&
                            trial.violation <= base.violation * 1.01 + 1e-7))
                        : (trial.merit < base.merit &&
                           (continuationDiagnostic || trial.violation < 1e-8));
                const bool take = filterOK && progress;
                Vector sourcePrediction;
                const bool sourceAvailable=m_reducedSourceRay && !restore && base.merit>1e-6 &&
                                           sourceModel(trialCoefficients,sourcePrediction);
                double sourceError=0.;
                if (sourceAvailable)
                {
                    auto difference=trial.r;
                    for (int j=0;j<dimension;++j) difference[j]-=sourcePrediction[j];
                    sourceError=dualNorm(difference);
                    emit("source_ray_prediction", {{"fraction",J::Number(fraction)},
                        {"predicted_norm",J::Number(dualNorm(sourcePrediction))},
                        {"actual_norm",J::Number(trial.merit)},
                        {"model_error",J::Number(sourceError)},
                        {"linear_model_error",J::Number(dualNorm(error))},
                        {"non_k_error",J::Number(trial.response.error)},
                        {"required",J::Number(innerTolerance)},
                        {"actual_displacement",J::Number(Length(trialDisplacement))},
                        {"radius",J::Number(radius)}});
                }
                if (!chordReuse && !decomposed && dualNorm(error) > .05 * std::max(1., base.merit))
                {
                    const auto baseD = component(state);
                    auto predictedD = baseD;
                    for (int i = 0; i < n; ++i)
                    {
                        const auto shiftedD = component(columns[i].shifted.response.state);
                        for (int j = 0; j < dimension; ++j)
                            predictedD[j] += trialCoefficients[i] *
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
                      {"proposal_origin",J::Quote(proposal.origin)},
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
                      {"second_order_corrections", J::Number(usedCorrections)},
                      {"displacement",
                       J::Number(Length(trialDisplacement))},
                      {"physical_k_displacement",
                       J::Number(kScale * Length(trialDisplacement))}});
                if (m_reducedVerify)
                {
                    const std::string file = "reduced-evaluated-" + stageTag + "-" +
                        std::to_string(iterations) + "-" + std::to_string(failedModels) +
                        "-p" + std::to_string(proposalIndex) + "-f" + J::Number(fraction) + ".json";
                    std::ofstream out(file);
                    out << J::Object({{"status",J::Quote("nonlinear_trial_only")},
                        {"state",J::Numbers(trial.response.state)},
                        {"history",J::Numbers(history)}, {"alpha",J::Number(alpha)},
                        {"time",J::Number(time)}, {"require_history_h0","false"},
                        {"inner_stress_model",J::Quote(InnerStressName())}}) << '\n';
                    out.close();
                    Require(out.good(),"failed to save evaluated reduced trial");
                    emit("evaluated_trial_packet", {{"path",J::Quote(file)},
                        {"proposal_origin",J::Quote(proposal.origin)},
                        {"fraction",J::Number(fraction)}, {"accepted",J::Boolean(take)}});
                }
                if (!take)
                    continue;
                const double oldViolation = base.violation;
                const double oldMerit = base.merit;
                if (m_reducedResponseReuseUpdates)
                {
                    try
                    {
                        chordColumns.clear();
                        // A secant model is temporary search information, never
                        // a partially equilibrated Y(k). Only converged actual
                        // endpoints supply changes, and each reused direction
                        // must pass new original-operator checks above.
                        if (!restore && trial.violation<1e-8 &&
                            trial.domain.profileBranches==base.domain.profileBranches &&
                            (!chordReuse || chordAge<m_reducedResponseReuseUpdates) &&
                            dualNorm(error)<.2*std::max(1.,base.merit))
                        {
                            RANSWallReduced::Matrix kBasis,rActions,stateActions,traceActions;
                            for (const auto &column:columns)
                            {
                                kBasis.push_back(column.k);rActions.push_back(column.r);
                                stateActions.push_back(column.stateAction);traceActions.push_back(column.traceAction);
                            }
                            Vector mapped(dimension),dr=trial.r,
                                   ds=trial.domain.effectiveState,dt=trial.domain.traceState;
                            for (int i=0;i<n;++i)
                                for (int j=0;j<dimension;++j) mapped[j]+=trialCoefficients[i]*columns[i].k[j];
                            for (int j=0;j<dimension;++j)
                            { mapped[j]-=trialDisplacement[j];dr[j]-=base.r[j]; }
                            for (int j=0;j<Size();++j)
                            { ds[j]-=base.domain.effectiveState[j];dt[j]-=base.domain.traceState[j]; }
                            const double spanError=Length(mapped);
                            const double minimumSignal=1e-8;
                            const bool updated=spanError<=1e-10*std::max(1.,Length(trialDisplacement)) &&
                                RANSWallReduced::UpdateSecantActions(rActions,kBasis,trialCoefficients,dr,minimumSignal) &&
                                RANSWallReduced::UpdateSecantActions(stateActions,kBasis,trialCoefficients,ds,minimumSignal) &&
                                RANSWallReduced::UpdateSecantActions(traceActions,kBasis,trialCoefficients,dt,minimumSignal);
                            double residualDefect=0.,stateDefect=0.,traceDefect=0.;
                            if (updated)
                            {
                                auto defect=[&](const RANSWallReduced::Matrix &actions,Vector change) {
                                    for (int i=0;i<n;++i)
                                        for (int j=0;j<change.size();++j) change[j]-=trialCoefficients[i]*actions[i][j];
                                    return change;
                                };
                                residualDefect=dualNorm(defect(rActions,dr))/std::max(1.,dualNorm(dr));
                                stateDefect=Length(defect(stateActions,ds))/std::max(1.,Length(ds));
                                traceDefect=Length(defect(traceActions,dt))/std::max(1.,Length(dt));
                                if (residualDefect>1e-8 || stateDefect>1e-8 || traceDefect>1e-8)
                                    throw std::runtime_error("inaccurate reduced response secant interpolation");
                                for (int i=0;i<n;++i)
                                {
                                    Column column=columns[i];
                                    column.r=std::move(rActions[i]);
                                    column.stateAction=std::move(stateActions[i]);
                                    column.traceAction=std::move(traceActions[i]);
                                    // k is an independent prescribed direction,
                                    // not an approximate part of the y response.
                                    for (int j=offset[K()];j<offset[Epsilon()];++j)
                                        column.stateAction[j]=column.traceAction[j]=0.;
                                    for (int j=0;j<dimension;++j)
                                        column.stateAction[rows[j]]=column.traceAction[rows[j]]=kScale*column.k[j];
                                    column.shifted=Evaluation{};
                                    column.observations.clear();column.observations.shrink_to_fit();
                                    chordColumns.push_back(std::move(column));
                                }
                                chordAge=chordReuse?chordAge+1:1;
                                if (!chordReuse) chordOrigin=iterations;
                                chordBranches=trial.domain.profileBranches;
                            }
                            emit("response_chord_secant",{{"updated",J::Boolean(updated)},
                                {"next_age",J::Number(chordAge)},
                                {"field_signal",J::Number(Length(trialDisplacement))},
                                {"minimum_field_signal",J::Number(minimumSignal)},
                                {"span_error",J::Number(spanError)},
                                {"residual_secant_defect",J::Number(residualDefect)},
                                {"state_secant_defect",J::Number(stateDefect)},
                                {"trace_secant_defect",J::Number(traceDefect)},
                                {"scope",J::Quote("secant interpolation defects only; not current derivative or nonlinear stationarity certificates")}});
                        }
                    }
                    catch (const std::runtime_error &error)
                    {
                        chordColumns.clear();chordAge=0;
                        emit("response_chord_discarded",{{"reason",J::Quote(error.what())},
                            {"scope",J::Quote("optional model update failed; retain the accepted actual trial and rebuild fresh responses")}});
                    }
                }
                previous = trialDisplacement;
                state               = trial.response.state;
                m_lastFeasibleState = state;
                base                = std::move(trial);
                filter.emplace_back(base.merit, base.violation);
                const double predictionError = dualNorm(error);
                const bool agreement = predictionError < .2 * std::max(1., base.merit)
                    && std::abs(base.violation - modelViolation) <
                       .2 * std::max(1., base.violation);
                if ((m_reducedSecondOrder || continuationDiagnostic) && !restore)
                {
                    const double predictedNorm=sourceAvailable ? dualNorm(sourcePrediction) : dualNorm(model);
                    const double predictedGain = oldMerit*oldMerit-predictedNorm*predictedNorm;
                    const double actualGain = oldMerit*oldMerit-base.merit*base.merit;
                    const double ratio = predictedGain > 0. ? actualGain/predictedGain : 0.;
                    if (ratio < .25) radius *= .5;
                    else if (ratio > .75 && Length(previous) > .8*radius)
                        radius = std::min(m_reducedMaximumRadius,1.8*radius);
                    emit("radius_agreement", {{"actual_to_predicted_gain", J::Number(ratio)},
                        {"prediction_model",J::Quote(sourceAvailable ? "nonlinear_destruction_surrogate" : "original_linear_reduced")},
                        {"source_model_evaluations",J::Number(sourceModelEvaluations)},
                        {"next_radius", J::Number(radius)}});
                }
                else
                    radius = std::min(m_reducedMaximumRadius,
                                      radius * (fraction == 1. && agreement ? 1.8 : .6));
                ++iterations;
                failedModels = 0;
                report(base, "accepted");
                directionRefreshes=0;
                if (m_reducedAdaptiveSpace && !restore && base.merit > 1e-6)
                {
                    const double relativeGain = (oldMerit-base.merit)/oldMerit;
                    if (fraction <= .25 || relativeGain < 1e-3)
                    {
                        const int oldBudget = basisBudget;
                        basisBudget = std::min(m_reducedMaximumBasis,
                            std::max(basisBudget+4,5*basisBudget/4));
                        const bool refresh = oldBudget == m_reducedMaximumBasis &&
                            (fraction <= .0625 || relativeGain < 1e-5);
                        if (refresh)
                        {
                            // Drop seed history only. Every governing row,
                            // actual qualification observation and physical
                            // history remains part of the next native solve.
                            activeWitnesses.clear();
                            previous.clear();
                        }
                        emit("response_space_adapt", {{"previous_budget",J::Number(oldBudget)},
                            {"next_budget",J::Number(basisBudget)},
                            {"seed_history_refreshed",J::Boolean(refresh)},
                            {"relative_original_residual_gain",J::Number(relativeGain)},
                            {"fraction",J::Number(fraction)}});
                    }
                }
                if (restore && base.violation > 0.)
                {
                    const double restored = (oldViolation - base.violation) /
                                            oldViolation;
                    qualificationStalls = restored < .01 ? qualificationStalls + 1 : 0;
                    if (qualificationStalls >= 2)
                    {
                        if (restorationWeight > 0.)
                            restorationWeight = std::max(1e-4, .1 * restorationWeight);
                        basisBudget = std::min(m_reducedMaximumBasis, basisBudget + 3);
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
            if (chordReuse)
            {
                chordColumns.clear();chordAge=0;
                emit("response_chord_discarded",{{"reason",J::Quote("no acceptable original-equation trial; rebuild fresh responses before changing radius")}});
                continue;
            }
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
    Require(continuationDiagnostic || audit.at("r2_qualified")[0] == 1.,
            "reduced algebraic result failed original qualification");
    return state;
}
} // namespace Nektar
